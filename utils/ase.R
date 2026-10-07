# Hospitalization-level ASE classification, adapted from pinned clifpy SQL (Apache-2.0).
# No Python runtime is required. Calendar windows use the project's UTC timestamps.
ase_required_columns <- list(
  hospitalization=c("hospitalization_id","patient_id","admission_dttm","discharge_dttm","discharge_category","age_at_admission"),
  patient=c("patient_id","death_dttm"),
  microbiology_culture=c("hospitalization_id","collect_dttm","fluid_category","method_category"),
  medication_admin_intermittent=c("hospitalization_id","admin_dttm","med_category","med_group","med_route_category","med_dose","mar_action_group"),
  medication_admin_continuous=c("hospitalization_id","admin_dttm","med_category","med_group","med_dose","mar_action_group"),
  adt=c("hospitalization_id","in_dttm","out_dttm","location_category"),
  labs=c("hospitalization_id","lab_category","lab_value","lab_value_numeric","lab_result_dttm","lab_order_dttm"),
  respiratory_support=c("hospitalization_id","recorded_dttm","device_category"),
  hospital_diagnosis=c("hospitalization_id","diagnosis_code")
)
ase_sql <- function(name, root="utils/ase/sql") paste(readLines(file.path(root,paste0(name,".sql")),warn=FALSE),collapse="\n")

ase_source_check <- function() {
  paths<-vapply(names(ase_required_columns),find_table_path,character(1),required=FALSE)
  missing<-names(paths)[is.na(paths)]
  if(length(missing))return(list(status="skipped_required_tables_unavailable",detail=paste(missing,collapse="; "),paths=paths))
  bad<-vapply(names(paths),function(n)paste(setdiff(ase_required_columns[[n]],janitor::make_clean_names(preflight_header(paths[[n]]))),collapse=","),character(1))
  if(any(nzchar(bad)))return(list(status="skipped_required_columns_unavailable",detail=paste(names(bad)[nzchar(bad)],bad[nzchar(bad)],sep=": ",collapse="; "),paths=paths))
  list(status="available",detail=NA_character_,paths=paths)
}

ase_connect <- function(temp_dir) {
  dir.create(temp_dir,recursive=TRUE,showWarnings=FALSE)
  con<-DBI::dbConnect(duckdb::duckdb(),dbdir=":memory:")
  DBI::dbExecute(con,"SET memory_limit='2GB'")
  DBI::dbExecute(con,"SET threads=4")
  DBI::dbExecute(con,paste("SET temp_directory=",DBI::dbQuoteString(con,normalizePath(temp_dir))))
  con
}

ase_register_sources <- function(con,paths,hospitalization_ids) {
  DBI::dbWriteTable(con,"cohort_ids",data.frame(hospitalization_id=as.character(unique(hospitalization_ids))))
  for(tbl in names(paths)) {
    columns<-ase_required_columns[[tbl]]
    path<-paths[[tbl]];ext<-tolower(tools::file_ext(path))
    if(ext %in% c("parquet","csv")) {
      reader<-if(ext=="parquet")"read_parquet" else "read_csv_auto"
      # Read CSV as text to preserve identifiers; timestamps/numerics are cast below.
      options<-if(ext=="csv")", header=true, all_varchar=true" else ""
      header<-preflight_header(path);raw<-header[match(columns,janitor::make_clean_names(header))]
      selection<-paste(paste(DBI::dbQuoteIdentifier(con,raw),"AS",DBI::dbQuoteIdentifier(con,columns)),collapse=",")
      query<-paste0("SELECT ",selection," FROM ",reader,"(",DBI::dbQuoteString(con,path),options,")")
      if(tbl!="patient")query<-paste0("SELECT t.* FROM (",query,") t JOIN cohort_ids c ON CAST(t.hospitalization_id AS VARCHAR)=c.hospitalization_id")
      DBI::dbExecute(con,paste0("CREATE VIEW raw_",tbl," AS ",query))
    } else {
      x<-preflight_columns(path,columns)
      key<-if(tbl=="patient")"patient_id" else "hospitalization_id"
      x[[key]]<-as.character(x[[key]])
      if(key=="hospitalization_id")x<-x[x[[key]] %in% hospitalization_ids,,drop=FALSE]
      DBI::dbWriteTable(con,paste0("raw_",tbl),x)
    }
  }
}

ase_standardize_sources <- function(con) {
  for(tbl in names(ase_required_columns)) {
    columns<-ase_required_columns[[tbl]]
    expressions<-vapply(columns,function(col){
      z<-as.character(DBI::dbQuoteIdentifier(con,col))
      if(grepl("_dttm$",col))paste0("TRY_CAST(",z," AS TIMESTAMP) AS ",z)
      else if(col %in% c("age_at_admission","med_dose","lab_value_numeric"))paste0("TRY_CAST(",z," AS DOUBLE) AS ",z)
      else if(grepl("_category$|_group$",col))paste0("lower(trim(CAST(",z," AS VARCHAR))) AS ",z)
      else paste0("CAST(",z," AS VARCHAR) AS ",z)
    },character(1))
    DBI::dbExecute(con,paste0("CREATE VIEW src_",tbl," AS SELECT ",paste(expressions,collapse=",")," FROM raw_",tbl))
  }
}

compute_ase_hospitalizations <- function(con,sql_root="utils/ase/sql") {
  exec<-function(sql)DBI::dbExecute(con,sql)
  # Only adults with complete hospitalization boundaries enter the two-group analysis.
  exec("CREATE TABLE hospitalizations AS SELECT h.* FROM src_hospitalization h JOIN cohort_ids c USING(hospitalization_id) WHERE h.age_at_admission>=18 AND h.admission_dttm IS NOT NULL AND h.discharge_dttm>h.admission_dttm")
  if(DBI::dbGetQuery(con,"SELECT count(*) n FROM hospitalizations")$n==0)return(list(classification=tibble::tibble(hospitalization_id=character(),ase_group=character(),presumed_infection=logical(),ase_without_lactate=logical(),ase_with_lactate=logical()),episodes=tibble::tibble()))
  if(DBI::dbGetQuery(con,"SELECT count(*) n FROM (SELECT hospitalization_id FROM hospitalizations GROUP BY hospitalization_id HAVING count(*)>1)")$n>0)stop("ASE: duplicate hospitalization IDs.")
  exec("CREATE VIEW patient AS SELECT p.* FROM src_patient p WHERE p.patient_id IN (SELECT patient_id FROM hospitalizations)")
  exec("CREATE TABLE blood_cultures AS SELECT m.hospitalization_id, m.collect_dttm AS culture_time, CAST(m.collect_dttm AS DATE) AS culture_day, row_number() OVER(PARTITION BY m.hospitalization_id ORDER BY m.collect_dttm) AS bc_id FROM (SELECT DISTINCT hospitalization_id,collect_dttm FROM src_microbiology_culture WHERE fluid_category='blood_buffy' AND method_category='culture' AND collect_dttm IS NOT NULL) m JOIN hospitalizations h USING(hospitalization_id) WHERE CAST(m.collect_dttm AS DATE)>=CAST(h.admission_dttm AS DATE)-2 AND m.collect_dttm<=h.discharge_dttm")
  exec("CREATE VIEW antibiotics AS SELECT m.*, CAST(admin_dttm AS DATE) AS med_admin_day, CASE WHEN med_route_category IN ('iv','im','intravenous','intramuscular') THEN 1 ELSE 0 END AS is_iv_im FROM src_medication_admin_intermittent m WHERE med_group='cms_sepsis_qualifying_antibiotics' AND mar_action_group='administered' AND med_dose>0 AND med_category IS NOT NULL")
  qad<-ase_sql("qad_query",sql_root)
  # Avoid losing admission-day treatment by comparing DATE with DATE. Include ED days.
  qad<-gsub("a.med_admin_day >= h.admission_dttm","a.med_admin_day >= CAST(h.admission_dttm AS DATE) - 2",qad,fixed=TRUE)
  qad<-gsub("a.med_admin_day <= h.discharge_dttm","a.med_admin_day <= CAST(h.discharge_dttm AS DATE)",qad,fixed=TRUE)
  exec(paste("CREATE TABLE qad_results AS",qad))
  exec(paste("CREATE TABLE final_qad AS",ase_sql("qad_censoring_query",sql_root)))
  exec("CREATE TABLE bc_episodes AS SELECT b.hospitalization_id,b.bc_id,b.culture_time AS blood_culture_dttm,b.culture_day AS blood_culture_day,COALESCE(q.meets_qad_with_censoring,0) AS meets_qad_with_censoring,q.anchor_meds_in_window,q.anchor_parenteral_meds_in_window,q.run_meds FROM blood_cultures b LEFT JOIN final_qad q USING(hospitalization_id,bc_id)")
  exec("CREATE VIEW labs AS SELECT * FROM src_labs")
  # N18.6 is the explicit ESRD diagnosis. Avoid upstream I27.2 (pulmonary hypertension).
  exec("CREATE TABLE esrd_patients AS SELECT DISTINCT hospitalization_id,1 AS has_esrd FROM src_hospital_diagnosis WHERE lower(replace(trim(diagnosis_code),'.',''))='n186'")
  exec(paste("CREATE TABLE lab_dysfunction AS",ase_sql("lab_dysfunction_query",sql_root)))
  exec("CREATE VIEW med_continuous AS SELECT * FROM src_medication_admin_continuous WHERE mar_action_group='administered' AND med_category IN ('norepinephrine','dopamine','epinephrine','phenylephrine','vasopressin')")
  exec("CREATE VIEW adt AS SELECT * FROM src_adt")
  exec("CREATE VIEW respiratory AS SELECT * FROM src_respiratory_support")
  exec("CREATE VIEW blood_cultures_temp AS SELECT hospitalization_id,bc_id,culture_time AS blood_culture_dttm FROM blood_cultures")
  exec(paste("CREATE TABLE vasopressor_df AS",ase_sql("vasopressor_query",sql_root)))
  exec(paste("CREATE TABLE imv_df AS",ase_sql("imv_query",sql_root)))
  exec(ase_sql("component_b_inputs",sql_root))
  episodes<-DBI::dbGetQuery(con,"SELECT *, presumed_infection=1 AND (vasopressor_dttm IS NOT NULL OR imv_dttm IS NOT NULL OR aki_dttm IS NOT NULL OR hyperbilirubinemia_dttm IS NOT NULL OR thrombocytopenia_dttm IS NOT NULL) AS ase_without_lactate, presumed_infection=1 AND (vasopressor_dttm IS NOT NULL OR imv_dttm IS NOT NULL OR aki_dttm IS NOT NULL OR hyperbilirubinemia_dttm IS NOT NULL OR thrombocytopenia_dttm IS NOT NULL OR lactate_dttm IS NOT NULL) AS ase_with_lactate FROM component_b_inputs")
  # Any qualifying episode assigns the hospitalization. RIT changes episode counts,
  # not this binary any-event classification; no repeated events are counted here.
  flags<-episodes %>% dplyr::group_by(hospitalization_id) %>% dplyr::summarise(presumed_infection=any(presumed_infection==1),ase_without_lactate=any(ase_without_lactate),ase_with_lactate=any(ase_with_lactate),.groups="drop")
  hospitals<-DBI::dbGetQuery(con,"SELECT hospitalization_id FROM hospitalizations")
  classification<-hospitals %>% dplyr::left_join(flags,by="hospitalization_id") %>% dplyr::mutate(dplyr::across(c(presumed_infection,ase_without_lactate,ase_with_lactate),~dplyr::coalesce(.x,FALSE)),ase_group=dplyr::if_else(ase_without_lactate,"ASE","Non-ASE"))
  list(classification=classification,episodes=episodes)
}

build_ase_group_aggregates <- function(data,classification,months) {
  if(anyDuplicated(classification$hospitalization_id))stop("ASE: hospitalization classification is not unique.")
  if(anyNA(classification$ase_group)||any(!classification$ase_group %in% c("ASE","Non-ASE")))stop("ASE: invalid classification; unknown must not be assigned non-ASE.")
  keys<-dplyr::select(classification,hospitalization_id,ase_group)
  stays<-dplyr::inner_join(data$icu_admissions,keys,by="hospitalization_id")
  events<-dplyr::inner_join(data$events,keys,by="hospitalization_id") %>% dplyr::mutate(calendar_month=lubridate::floor_date(collect_dttm,"month"))
  rows<-dplyr::inner_join(data$rows,keys,by="hospitalization_id") %>% dplyr::mutate(calendar_month=lubridate::floor_date(collect_dttm,"month"))
  den<-dplyr::bind_rows(lapply(c("ASE","Non-ASE"),function(g)dplyr::mutate(monthly_icu_denominators(dplyr::filter(stays,ase_group==g),months),ase_group=g)))
  monthly<-events %>% dplyr::group_by(ase_group,calendar_month) %>% dplyr::summarise(n_culture_events=dplyr::n(),n_positive_culture_events=sum(any_positive_culture),.groups="drop")
  monthly<-den %>% dplyr::left_join(monthly,by=c("ase_group","calendar_month")) %>% dplyr::mutate(dplyr::across(c(n_culture_events,n_positive_culture_events),~dplyr::coalesce(.x,0L)),culture_events_per_100_icu_days=dplyr::if_else(n_icu_days>0,100*n_culture_events/n_icu_days,NA_real_),positive_events_per_100_icu_days=dplyr::if_else(n_icu_days>0,100*n_positive_culture_events/n_icu_days,NA_real_),culture_events_per_100_icu_admissions=dplyr::if_else(n_icu_admissions>0,100*n_culture_events/n_icu_admissions,NA_real_),positivity=dplyr::if_else(n_culture_events>0,n_positive_culture_events/n_culture_events,NA_real_))
  # Use pooled site source activity, so a subgroup with cultures=0 has a real zero
  # when cultures continue elsewhere at the site. Never condition exposure on culturing.
  source<-data$events %>% dplyr::count(calendar_month=lubridate::floor_date(collect_dttm,"month"),name="n_observed_culture_events")
  monthly<-monthly %>% dplyr::left_join(source,by="calendar_month") %>% dplyr::mutate(n_observed_culture_events=dplyr::coalesce(n_observed_culture_events,0L))
  organisms<-rows %>% dplyr::filter(positive_culture,!is.na(organism_category),!organism_category %in% c("na","unknown","missing")) %>% dplyr::distinct(ase_group,calendar_month,culture_event_id,organism_category) %>% dplyr::count(ase_group,calendar_month,organism_category,name="n_detection_events")
  organisms<-tidyr::crossing(dplyr::select(monthly,ase_group,calendar_month,n_icu_days,n_icu_admissions,n_observed_culture_events),organism_category=sort(unique(organisms$organism_category))) %>% dplyr::left_join(organisms,by=c("ase_group","calendar_month","organism_category")) %>% dplyr::mutate(n_detection_events=dplyr::coalesce(n_detection_events,0L),detections_per_100_icu_days=dplyr::if_else(n_icu_days>0,100*n_detection_events/n_icu_days,NA_real_))
  specimen<-events %>% dplyr::mutate(fluid_category=dplyr::coalesce(fluid_category,"missing")) %>% dplyr::group_by(ase_group,calendar_month,fluid_category) %>% dplyr::summarise(n_culture_events=dplyr::n(),n_positive_culture_events=sum(any_positive_culture),.groups="drop")
  specimen<-tidyr::crossing(den,fluid_category=sort(unique(specimen$fluid_category))) %>% dplyr::left_join(specimen,by=c("ase_group","calendar_month","fluid_category")) %>% dplyr::mutate(dplyr::across(c(n_culture_events,n_positive_culture_events),~dplyr::coalesce(.x,0L)),cultures_per_100_icu_days=dplyr::if_else(n_icu_days>0,100*n_culture_events/n_icu_days,NA_real_),positivity=dplyr::if_else(n_culture_events>0,n_positive_culture_events/n_culture_events,NA_real_))
  first<-events %>% dplyr::group_by(icu_admission_id) %>% dplyr::summarise(first_culture_dttm=min(collect_dttm),.groups="drop")
  admission<-stays %>% dplyr::filter(icu_in_dttm>=icu_in_dttm_clipped) %>% dplyr::left_join(first,by="icu_admission_id") %>% dplyr::mutate(hours_to_first_culture=as.numeric(difftime(first_culture_dttm,icu_in_dttm,units="hours")),cultured=!is.na(first_culture_dttm))
  timing<-admission %>% dplyr::group_by(ase_group) %>% dplyr::summarise(n_icu_admissions=dplyr::n(),n_with_icu_culture=sum(cultured),percent_with_icu_culture=100*mean(cultured),median_hours_to_first_culture=if(any(cultured))median(hours_to_first_culture[cultured]) else NA_real_,p25_hours=if(any(cultured))unname(quantile(hours_to_first_culture[cultured],.25)) else NA_real_,p75_hours=if(any(cultured))unname(quantile(hours_to_first_culture[cultured],.75)) else NA_real_,.groups="drop")
  cohort<-stays %>% dplyr::group_by(ase_group) %>% dplyr::summarise(n_patients=dplyr::n_distinct(patient_id),n_hospitalizations=dplyr::n_distinct(hospitalization_id),n_icu_stays=dplyr::n(),n_icu_days=sum(icu_los_days),.groups="drop")
  admission_monthly<-admission %>% group_by(ase_group,calendar_month=icu_admission_month) %>% summarise(n_cultured_icu_admissions=sum(cultured),.groups="drop")
  admission_monthly<-den %>% left_join(admission_monthly,by=c("ase_group","calendar_month")) %>% mutate(n_cultured_icu_admissions=coalesce(n_cultured_icu_admissions,0L),fraction_admissions_cultured=if_else(n_icu_admissions>0,n_cultured_icu_admissions/n_icu_admissions,NA_real_))
  bins<-admission %>% mutate(timing_bin=case_when(!cultured~"No ICU culture",hours_to_first_culture<=6~"0-6 hours",hours_to_first_culture<=24~"6-24 hours",hours_to_first_culture<=48~"24-48 hours",hours_to_first_culture<=168~"2-7 days",TRUE~">7 days")) %>% count(ase_group,timing_bin,name="n_icu_admissions") %>% group_by(ase_group) %>% mutate(percent_admissions=100*n_icu_admissions/sum(n_icu_admissions)) %>% ungroup()
  list(stays=stays,rows=rows,events=events,monthly=monthly,organisms=organisms,specimen=specimen,timing=timing,cohort=cohort,admission_monthly=admission_monthly,timing_bins=bins)
}
