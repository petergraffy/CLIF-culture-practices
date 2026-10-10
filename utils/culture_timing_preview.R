# Scoped culture loading and ICU location for the antibiotic timing preview.
read_timing_cultures <- function(con, ids) {
  path<-find_table_path("microbiology_culture")
  if(tolower(tools::file_ext(path)) %in% c("parquet","csv")) {
    DBI::dbWriteTable(con,"culture_timing_ids",data.frame(hospitalization_id=as.character(ids)),overwrite=TRUE)
    reader<-if(tolower(tools::file_ext(path))=="parquet")"read_parquet" else "read_csv_auto"
    options<-if(reader=="read_csv_auto")",header=true,all_varchar=true" else ""
    header<-preflight_header(path);clean<-janitor::make_clean_names(header)
    cols<-intersect(c("patient_id","hospitalization_id","organism_id","organism_name","organism_category","organism_group","order_dttm","collect_dttm","result_dttm","fluid_name","fluid_category","method_name","method_category"),clean)
    selected<-vapply(cols,function(n)paste0(if(n %in% c("patient_id","hospitalization_id","organism_id"))"CAST(" else "",DBI::dbQuoteIdentifier(con,header[match(n,clean)]),if(n %in% c("patient_id","hospitalization_id","organism_id"))" AS VARCHAR)" else ""," AS ",DBI::dbQuoteIdentifier(con,n)),character(1))
    query<-paste0("SELECT t.* FROM (SELECT ",paste(selected,collapse=",")," FROM ",reader,"(",DBI::dbQuoteString(con,path),options,")) t JOIN culture_timing_ids i USING(hospitalization_id)")
    DBI::dbGetQuery(con,query)
  } else read_tbl("microbiology_culture") %>% mutate(hospitalization_id=as.character(hospitalization_id)) %>% filter(hospitalization_id %in% ids)
}

anchor_antibiotic_phase <- function(markers, stays, start, end) {
  first<-markers %>% arrange(marker_time,marker_type) %>% group_by(hospitalization_id) %>% slice_head(n=1) %>% ungroup()
  limits<-stays %>% group_by(hospitalization_id) %>% summarise(first_icu_in=min(icu_in_dttm),last_icu_out=max(icu_out_dttm),.groups="drop")
  matched<-first %>% inner_join(stays %>% select(hospitalization_id,icu_in_dttm,icu_out_dttm),by="hospitalization_id",relationship="many-to-many") %>% filter(marker_time>=icu_in_dttm,marker_time<icu_out_dttm) %>% select(hospitalization_id,icu_in_dttm)
  if(anyDuplicated(matched$hospitalization_id))stop("Support marker matches multiple ICU intervals")
  first %>% left_join(limits,by="hospitalization_id") %>% left_join(matched,by="hospitalization_id") %>% mutate(hours_into_icu=as.numeric(difftime(marker_time,icu_in_dttm,units="hours")),anchor_phase=case_when(marker_time<first_icu_in~"Before first ICU entry",!is.na(icu_in_dttm)&hours_into_icu<24~"First 24 ICU hours",!is.na(icu_in_dttm)~"After 24 ICU hours",TRUE~"Between/after ICU stays"),within_study=marker_time>=start&marker_time<end)
}
