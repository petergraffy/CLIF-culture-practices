# Startup checks use aggregate metadata only; no patient identifiers are exported.
preflight_packages <- function(file_type = NULL, check_versions = TRUE) {
  packages <- c("jsonlite","dplyr","tidyr","readr","lubridate","ggplot2","glue","janitor","stringr","forcats","scales","purrr","mgcv","DBI","duckdb")
  packages <- unique(c(packages,if(identical(file_type,"parquet"))"arrow",if(identical(file_type,"fst"))"fst"))
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Preflight: missing packages: ",paste(missing,collapse=", "),". From the repository root, run Rscript -e 'renv::restore(prompt = FALSE)' and retry.",call.=FALSE)
  if(check_versions) {
    lock <- jsonlite::fromJSON("renv.lock")
    mismatched <- packages[vapply(packages,function(p) is.null(lock$Packages[[p]]) || packageVersion(p)!=package_version(lock$Packages[[p]]$Version),logical(1))]
    if(length(mismatched))stop("Preflight: package versions differ from renv.lock: ",paste(mismatched,collapse=", "),". Run renv::restore(prompt = FALSE) and restart R.",call.=FALSE)
  }
  invisible(packages)
}

preflight_header <- function(path) {
  switch(tolower(tools::file_ext(path)),
    csv=names(readr::read_csv(path,n_max=0,show_col_types=FALSE)),
    parquet=arrow::open_dataset(path,format="parquet")$schema$names,
    fst=fst::metadata_fst(path)$columnNames,
    stop("Preflight: unsupported table format: ",path,call.=FALSE))
}

preflight_columns <- function(path, columns) {
  # Preserve the same cleaned column names as read_any(), while reading only needed columns.
  header <- preflight_header(path)
  raw <- header[match(columns,janitor::make_clean_names(header))]
  x <- switch(tolower(tools::file_ext(path)),
    csv=readr::read_csv(path,col_types=do.call(readr::cols_only,setNames(rep(list(readr::col_character()),length(raw)),raw)),show_col_types=FALSE),
    parquet=arrow::read_parquet(path,col_select=dplyr::all_of(raw)),
    fst=fst::read_fst(path,columns=raw))
  names(x)<-janitor::make_clean_names(names(x))
  x
}

preflight_schema <- function(table,header,ast_available=FALSE) {
  required <- switch(table,
    hospitalization=c("patient_id","hospitalization_id","admission_dttm","discharge_dttm"),
    adt=c("hospitalization_id","in_dttm","out_dttm","location_category"),
    microbiology_culture=c("patient_id","hospitalization_id","collect_dttm","method_category","fluid_name","fluid_category","organism_name","organism_category","organism_group",if(ast_available)"organism_id"),
    microbiology_susceptibility=c("organism_id","antimicrobial_category","susceptibility_category"))
  missing <- setdiff(required,janitor::make_clean_names(header))
  if(length(missing))stop("Preflight: ",table," is missing required columns: ",paste(missing,collapse=", "),". Check the CLIF export schema.",call.=FALSE)
  invisible(TRUE)
}

preflight_date_summary <- function(table,x,start,end) {
  if(table=="microbiology_culture") {
    z<-safe_ts(x$collect_dttm);keep<-clean_micro_label(x$method_category)=="culture";keep[is.na(keep)]<-FALSE
    z<-z[keep];in_window<-!is.na(z)&z>=start&z<end;dates<-z
  } else {
    if(table=="adt") {keep<-clean_micro_label(x$location_category)=="icu";keep[is.na(keep)]<-FALSE;x<-x[keep,,drop=FALSE]}
    a<-safe_ts(x[[if(table=="adt")"in_dttm" else "admission_dttm"]]);b<-safe_ts(x[[if(table=="adt")"out_dttm" else "discharge_dttm"]])
    in_window<-!is.na(a)&a<end&(is.na(b)|b>start);dates<-c(a,b)
  }
  dates<-dates[!is.na(dates)]
  if(!length(dates) || !any(in_window))stop("Preflight: no usable ",table,if(table=="adt")" ICU" else ""," records overlap the shared study window. Check the source tables.",call.=FALSE)
  tibble::tibble(table=table,available_start_date=as.Date(min(dates)),available_end_date=as.Date(max(dates)),n_source_rows_in_study_window=sum(in_window),covers_start_month=lubridate::floor_date(min(dates),"month")<=lubridate::floor_date(start,"month"),covers_end_month=lubridate::floor_date(max(dates),"month")>=lubridate::floor_date(end-lubridate::seconds(1),"month"))
}

run_site_preflight <- function() {
  if(!dir.exists(clif_repo_path))stop("Preflight: repo directory does not exist. Set repo to the local repository's absolute path.",call.=FALSE)
  if(file.access(clif_repo_path,2)!=0)stop("Preflight: repo directory is not writable.",call.=FALSE)
  if(!dir.exists(clif_tables_path))stop("Preflight: tables_path directory does not exist.",call.=FALSE)
  if(!clif_file_type %in% c("csv","parquet","fst","auto"))stop("Preflight: file_type must be csv, parquet, fst, or auto.",call.=FALSE)
  if(!grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$",clif_site_name))stop("Preflight: site_name must use letters, numbers, underscores, or hyphens.",call.=FALSE)
  start<-safe_ts(study_settings$study_start_date);end<-safe_ts(study_settings$study_end_date)+lubridate::days(1)
  if(is.na(start)||is.na(end)||end<=start)stop("Preflight: invalid shared study dates in utils/study_settings.R.",call.=FALSE)
  paths<-vapply(c("hospitalization","adt","microbiology_culture"),find_table_path,character(1))
  ast<-find_table_path("microbiology_susceptibility",required=FALSE)
  if(!is.na(ast))paths<-c(paths,microbiology_susceptibility=ast)
  for(ext in unique(tolower(tools::file_ext(paths))))preflight_packages(ext)
  # Check every schema before scanning date columns, including optional AST when present.
  for(table in names(paths))preflight_schema(table,preflight_header(paths[[table]]),!is.na(ast))
  summaries<-lapply(names(paths)[names(paths)!="microbiology_susceptibility"],function(table){
    columns<-switch(table,hospitalization=c("admission_dttm","discharge_dttm"),adt=c("in_dttm","out_dttm","location_category"),microbiology_culture=c("collect_dttm","method_category"))
    preflight_date_summary(table,preflight_columns(paths[[table]],columns),start,end)
  })
  summary<-dplyr::mutate(dplyr::bind_rows(summaries),site_name=clif_site_name)
  message("Preflight passed. AST table: ",if(is.na(ast))"unavailable (analysis will skip)" else "available")
  print(summary)
  if(any(!summary$covers_start_month|!summary$covers_end_month))message("Source records do not span every boundary month of the shared window; report the available dates to the study team. The pipeline will use observed records.")
  list(summary=summary,ast_available=!is.na(ast))
}
