suppressPackageStartupMessages({library(dplyr);library(tidyr);library(lubridate)})
source("utils/culture_core.R");source("utils/ase_onset_preview.R")
t0<-as.POSIXct("2020-01-03 00:00:00",tz="UTC")
a<-tibble(hospitalization_id=c("A","B"),timing_anchor="Earliest qualifying organ dysfunction",event_time=t0,first_icu_in=t0+c(1,-1)*3600,icu_in_dttm=as.POSIXct(c(NA,as.numeric(t0)-3600),origin="1970-01-01",tz="UTC"))
h<-tibble(hospitalization_id=c("A","B"),admission_dttm=t0+c(-1,-49)*3600,discharge_dttm=t0+c(.5,73)*3600)
e<-tibble(hospitalization_id=c("A","B","B","B","B"),collect_dttm=t0+c(0,-48,0,6,72)*3600,fluid_category=c("blood_buffy","blood_buffy","blood_buffy","genito_urinary_tract","blood_buffy"),result_status=c("Positive","Positive","Positive","Negative/no growth","Positive"))
def<-tibble(hospitalization_id=c("A","B"),collect_dttm=t0)
x<-build_ase_onset_preview(e,a,h,def,t0-100*3600,t0+100*3600)
# Partial hospitalization contributes only half an observed hour in the [0,1) bin.
stopifnot(filter(x$exposure,onset_stratum=="At/before first ICU entry",bin_start_hour==0)$observed_patient_hours==.5)
# Left window boundary included, right boundary excluded; bin-edge collection at 6 is in [6,7).
y<-filter(x$yield,specimen=="All cultures",sensitivity=="All collections")
stopifnot(sum(y$n_culture_events)==4,sum(y$n_positive)==3,!anyNA(y$n_positive))
stopifnot(filter(y,onset_stratum=="During ICU stay",bin_start_hour==6)$n_culture_events==1)
y2<-filter(x$yield,specimen=="All cultures",sensitivity=="Exclude linked ASE blood collection")
stopifnot(sum(y2$n_culture_events)==2,sum(y2$n_positive)==1)
# Zero-collection bins have unavailable positivity, not false zero yield.
stopifnot(all(is.na(y$percent_positive[y$n_culture_events==0])))
cat("Onset preview boundary, partial-observation, positivity, and defining-collection checks passed.\n")
