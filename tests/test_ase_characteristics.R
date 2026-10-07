suppressPackageStartupMessages({library(dplyr);library(tidyr);library(lubridate)})
source("utils/culture_core.R");source("utils/ase_characteristics.R")
h <- tibble(hospitalization_id=c("a","b","c"),patient_id=c("p","p","q"),
  admission_dttm=safe_ts("2020-01-01"),discharge_dttm=safe_ts(c("2020-01-11","2020-01-21","2020-01-31")),
  age_at_admission=c(40,50,NA),admission_type_category=c("emergency","emergency",NA),discharge_category=c("home","expired","home"))
stays <- tibble(hospitalization_id=c("a","a","b","c"),icu_los_days=c(2,3,4,1))
events <- tibble(hospitalization_id=c("a","a"),culture_event_id=c("1","2"),any_positive_culture=c(TRUE,FALSE))
p <- tibble(patient_id=c("p","q"),sex_category=c("female","unknown"),race_category=c("white",NA))
x <- build_ase_characteristics(list(hospitalization=h,icu_admissions=stays,events=events),
  tibble(hospitalization_id=c("a","b","c"),ase_group=c("ASE","Non-ASE","Non-ASE")),p)
get <- function(field,group="Non-ASE",level="")filter(x$long,characteristic==field,ase_group==group,.data$level==.env$level)
stopifnot(get("age_at_admission")$median==50,get("age_at_admission")$n_observed==1,get("age_at_admission")$n_missing==1)
stopifnot(get("icu_days","ASE")$median==5,get("icu_stays","ASE")$median==2,get("hospitalizations","ASE")$n==1)
stopifnot(get("sex_category",level="female")$percent==50,get("sex_category",level="Missing / unknown")$percent==50)
stopifnot(get("ethnicity_category",level="Missing / unknown")$availability=="source_column_unavailable")
stopifnot(get("any_icu_culture",level="no")$n==2,get("any_positive_icu_culture","ASE",level="yes")$n==1)
stopifnot(x$qc$n[x$qc$metric=="patients_in_both_subgroups"]==1,x$qc$n[x$qc$metric=="unique_patients_overall"]==2)
stopifnot(!any(c("patient_id","hospitalization_id","admission_dttm") %in% names(x$long)))
empty <- build_ase_characteristics(list(hospitalization=h,icu_admissions=stays,events=events),
  tibble(hospitalization_id=c("a","b","c"),ase_group="ASE"),p)
stopifnot(all(filter(empty$long,ase_group=="Non-ASE")$availability=="empty_subgroup"))
conflict <- bind_rows(p,mutate(p[1,],sex_category="male"))
stopifnot(inherits(try(build_ase_characteristics(list(hospitalization=h,icu_admissions=stays,events=events),
  tibble(hospitalization_id=c("a","b","c"),ase_group="ASE"),conflict),silent=TRUE),"try-error"))
clinical <- tibble(hospitalization_id=c("a","b","c"),recorded_imv=c("yes","no","no"),recorded_vasopressor=c("yes","yes","no"))
supported <- build_ase_characteristics(list(hospitalization=h,icu_admissions=stays,events=events),
  tibble(hospitalization_id=c("a","b","c"),ase_group=c("ASE","Non-ASE","Non-ASE")),p,clinical)
stopifnot(filter(supported$long,characteristic=="recorded_vasopressor",ase_group=="Non-ASE",level=="yes")$percent==50)
stopifnot(filter(supported$long,characteristic=="in_hospital_death",ase_group=="Non-ASE",level=="yes")$n==1)
cat("ASE characteristics tests passed: hospitalization weighting, missingness, denominators, repeat patients, empty groups, duplicate protection.\n")
