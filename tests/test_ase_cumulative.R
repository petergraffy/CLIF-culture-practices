suppressPackageStartupMessages({library(dplyr);library(tidyr);library(lubridate)})
source("utils/culture_core.R");source("utils/ase.R");source("utils/ase_cumulative.R")
t0<-safe_ts("2020-01-01")
stays<-tibble(icu_admission_id=c("a","u","b","c","carry"),ase_group=ase_group_levels[c(1,1,2,3,3)],
  icu_in_dttm=t0+hours(c(0,0,0,0,-24)),icu_in_dttm_clipped=t0,icu_los_days=c(2,.25,9,9,2))
events<-tibble(icu_admission_id=c("a","a","a","b","b","carry"),culture_event_id=as.character(1:6),
  collect_dttm=t0+3600*c(0,24,24.01,168,169,1),any_positive_culture=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),fluid_category="blood_buffy") %>% left_join(select(stays,icu_admission_id,ase_group,icu_in_dttm),by="icu_admission_id")
rows<-filter(events,any_positive_culture) %>% mutate(positive_culture=TRUE,organism_category="escherichia_coli")
rows<-bind_rows(rows,rows[1,]) # Duplicate isolate row and repeat cultures cannot inflate first detections.
z<-list(stays=stays,events=events,rows=rows)
x<-build_ase_cumulative(z)
first<-function(day,outcome="First culture",g=ase_group_levels[1])filter(x$first,ase_group==g,icu_day==day,.data$outcome==.env$outcome)
stopifnot(first(0)$n_admissions_with_event_by_day==1,first(0)$cumulative_percent==50,
  first(1,"First positive culture")$cumulative_percent==50,first(0,"First positive culture")$cumulative_percent==0)
stopifnot(filter(x$organisms,ase_group==ase_group_levels[1],icu_hour==24)$n_admissions_with_organism_by_hour==1)
stopifnot(filter(x$cultures,ase_group==ase_group_levels[1],icu_hour==24)$n_culture_events_by_hour==2,
  filter(x$cultures,ase_group==ase_group_levels[1],icu_hour==25)$cumulative_events_per_100_admissions==150)
stopifnot(filter(x$organisms,ase_group==ase_group_levels[2],icu_hour==167)$cumulative_percent==0,
  filter(x$organisms,ase_group==ase_group_levels[2],icu_hour==168)$cumulative_percent==100)
stopifnot(all(filter(x$organisms,ase_group=="ASE")$cumulative_percent==0),sum(x$qc$n_icu_admissions)==4,
  sum(x$qc$n_carry_in_stays_excluded)==1)
stopifnot(filter(x$observation,ase_group==ase_group_levels[1],icu_hour==6)$n_icu_admissions_under_observation==1,
  filter(x$observation,ase_group==ase_group_levels[1],icu_hour==48)$n_icu_admissions_under_observation==0)
stopifnot(all(x$first$n_admissions_with_event_by_day<=x$first$n_icu_admissions),
  all(x$organisms$n_admissions_with_organism_by_hour<=x$organisms$n_icu_admissions))
empty<-z;empty$stays<-filter(stays,ase_group!=ase_group_levels[2]);empty$events<-filter(events,ase_group!=ase_group_levels[2]);empty$rows<-filter(rows,ase_group!=ase_group_levels[2])
e<-build_ase_cumulative(empty);stopifnot(all(is.na(filter(e$first,ase_group==ase_group_levels[2])$cumulative_percent)))
negative<-z;negative$rows<-rows[0,];n<-build_ase_cumulative(negative);stopifnot(nrow(n$organisms)==0)
no_events<-z;no_events$events<-events[0,];no_events$rows<-rows[0,]
none<-build_ase_cumulative(no_events);stopifnot(nrow(none$cultures)==0,nrow(none$organisms)==0,all(none$first$cumulative_percent==0))
stopifnot(inherits(try(build_ase_cumulative(z,max_hour=0),silent=TRUE),"try-error"))
cat("Cumulative tests passed: time-zero and endpoint events, repeated detections, uncultured denominators, carry-in exclusions, observation counts, empty groups and no positive organisms.\n")
