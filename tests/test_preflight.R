suppressPackageStartupMessages({library(testthat);library(dplyr);library(lubridate)})
source("utils/preflight.R");source("utils/culture_core.R")
test_that("schema failures identify the table and fields before analysis", {
 expect_error(preflight_schema("hospitalization",c("patient_id","hospitalization_id")),"admission_dttm, discharge_dttm")
 expect_error(preflight_schema("microbiology_susceptibility",c("organism_id","antimicrobial_category")),"susceptibility_category")
 header<-c("patient_id","hospitalization_id","collect_dttm","method_category","fluid_name","fluid_category","organism_name","organism_category","organism_group")
 expect_silent(preflight_schema("microbiology_culture",header,FALSE))
 expect_error(preflight_schema("microbiology_culture",header,TRUE),"organism_id")
})
start<-safe_ts("2018-01-01");end<-safe_ts("2025-01-01")
test_that("date checks require eligible source rows and retain carry-in ICU stays", {
 culture<-tibble(collect_dttm=c("2019-01-01","2025-01-01"),method_category="culture")
 x<-preflight_date_summary("microbiology_culture",culture,start,end)
 expect_equal(x$n_source_rows_in_study_window,1L);expect_equal(x$available_end_date,as.Date("2025-01-01"))
 expect_error(preflight_date_summary("microbiology_culture",mutate(culture,method_category="pcr"),start,end),"no usable")
 expect_error(preflight_date_summary("microbiology_culture",mutate(culture,collect_dttm="2025-02-01"),start,end),"shared study window")
 adt<-tibble(in_dttm="2017-12-30",out_dttm="2018-01-03",location_category="icu")
 expect_equal(preflight_date_summary("adt",adt,start,end)$n_source_rows_in_study_window,1L)
 expect_error(preflight_date_summary("adt",mutate(adt,location_category="ward"),start,end),"ICU records")
})
test_that("projected schema and date reads work for all supported formats", {
 dir<-tempfile("preflight_formats_");dir.create(dir)
 x<-tibble(collect_dttm=c("2019-01-01","2020-01-01"),method_category="culture",patient_id=c("PRIVATE_A","PRIVATE_B"))
 readr::write_csv(x,file.path(dir,"table.csv"));arrow::write_parquet(x,file.path(dir,"table.parquet"));fst::write_fst(as.data.frame(x),file.path(dir,"table.fst"))
 for(ext in c("csv","parquet","fst")) {
   path<-file.path(dir,paste0("table.",ext))
   expect_setequal(preflight_header(path),names(x))
   y<-preflight_columns(path,c("collect_dttm","method_category"))
   expect_setequal(names(y),c("collect_dttm","method_category"))
   expect_equal(preflight_date_summary("microbiology_culture",y,start,end)$n_source_rows_in_study_window,2L)
   preflight_packages(ext)
 }
})

test_that("missing dependencies produce an actionable early error", {
 env<-new.env(parent=globalenv());sys.source("utils/preflight.R",envir=env)
 env$requireNamespace<-function(package,quietly=TRUE)package!="mgcv"
 expect_error(env$preflight_packages(check_versions=FALSE),"missing packages: mgcv.*restore")
 env$requireNamespace<-function(package,quietly=TRUE)package!="fst"
 expect_silent(env$preflight_packages("csv",check_versions=FALSE))
 expect_error(env$preflight_packages("fst",check_versions=FALSE),"missing packages: fst.*restore")
})
