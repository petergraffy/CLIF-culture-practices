suppressPackageStartupMessages({library(dplyr);library(tidyr);library(lubridate);library(testthat)})
source("utils/culture_core.R");source("utils/susceptibility.R");source("utils/trends.R")
ts <- function(x) as.POSIXct(x, tz = "UTC")
h <- tibble(patient_id="p",hospitalization_id="h",admission_dttm=ts("2020-01-01"),discharge_dttm=ts("2020-02-05"))
adt <- tibble(hospitalization_id="h",location_category="icu",in_dttm=ts(c("2020-01-30", "2020-01-31", "2020-02-01", "2020-02-04")),out_dttm=ts(c("2020-02-01", "2020-02-02", "2020-02-03", "2020-02-05")))
micro <- tibble(patient_id="p", hospitalization_id="h", organism_id=c("a","a","b","c","d"), order_dttm=ts("2020-01-31"),collect_dttm=ts(c("2020-02-01","2020-02-01","2020-02-01","2020-02-03","2020-02-04")),result_dttm=ts(c("2020-02-01","2020-02-02","2020-02-02","2020-02-03","2020-02-04")),fluid_name="blood",fluid_category="blood_buffy",method_name="culture",method_category="culture",organism_name=c("S. aureus","S. aureus","E. coli","S. aureus",NA),organism_category=c("staphylococcus_aureus","staphylococcus_aureus","escherichia_coli","staphylococcus_aureus",NA),organism_group=c("staphylococcus","staphylococcus","escherichia","staphylococcus",NA))
test_that("overlapping and contiguous ADT rows merge and stop boundary is exclusive", {
 x<-build_culture_data(h,adt,micro);expect_equal(nrow(x$icu_admissions),2L);expect_equal(nrow(x$events),2L);expect_equal(nrow(x$rows),3L);expect_false("c" %in% x$rows$organism_id);expect_equal(length(unique(x$rows$culture_event_id)),2L)
 expect_equal(x$rows$result_dttm[x$rows$organism_id=="a"],ts("2020-02-02"))
})
test_that("unknown, mixed flora, and no growth are not called positive", {
 expect_equal(classify_culture_result(c(NA,"no_growth","mixed_flora","staphylococcus"),c(NA,"no_growth","mixed_flora","staphylococcus_aureus"),c(NA,"no growth","mixed flora","S. aureus")),culture_result_levels[c(4,1,3,2)])
 expect_equal(collapse_result_status(c("Negative/no growth","Indeterminate")),"Indeterminate")
})
test_that("carry-in stays contribute days but not new admissions", {
 x<-build_culture_data(h,adt,micro,start=ts("2020-02-01"),end=ts("2020-02-05"));den<-monthly_icu_denominators(x$icu_admissions,ts("2020-02-01"));expect_equal(den$n_icu_admissions,1L);expect_equal(den$n_icu_days,3)
})
test_that("specimen overrides are specific to source label and dates", {
 rows<-tibble(fluid_name="bal",fluid_category="other_unspecified",collect_dttm=ts(c("2020-01-01","2021-01-01")))
 rules<-tibble(fluid_name="bal",original_category="other_unspecified",start_date="2021-01-01",end_date="2021-12-31",replacement_category="respiratory_tract_lower",reason="source specimen verified")
 expect_equal(apply_specimen_overrides(rows,rules)$fluid_category,c("other_unspecified","respiratory_tract_lower"))
})
test_that("mCIDE categories drive AST; duplicate and conflicting tests do not inflate counts", {
 x<-build_culture_data(h,adt,micro);den<-monthly_icu_denominators(x$icu_admissions,ts(c("2020-01-01","2020-02-01")))
 ast<-tibble(organism_id=c("a","a","a","b","b","a"),antimicrobial_category=c(rep("oxacillin",3),"ceftriaxone","meropenem","vancomycin"),susceptibility_category=c("susceptible","susceptible","non_susceptible","non_susceptible",NA,NA),susceptibility_name="resistant")
 result<-build_susceptibility_analysis(x$rows,ast,den)
 ox<-filter(result$monthly,calendar_month==ts("2020-02-01"),antimicrobial_category=="oxacillin",specimen_stratum=="Overall")
 expect_equal(ox$n_indeterminate,1L);expect_equal(ox$n_interpretable,0L);expect_true(is.na(ox$non_susceptible_fraction));expect_equal(ox$n_conflicting,1L)
 mero<-filter(result$monthly,calendar_month==ts("2020-02-01"),antimicrobial_category=="meropenem",specimen_stratum=="Overall");expect_equal(mero$n_unavailable_reported,1L);expect_equal(mero$n_non_susceptible,0L)
 expect_error(normalize_susceptibility(tibble(organism_id="a")),"mCIDE fields")
 duplicate<-bind_rows(x$rows,mutate(x$rows[1,],culture_event_id=999L));expect_error(build_susceptibility_analysis(duplicate,ast,den),"multiple ICU")
})
test_that("NB and tested-fraction GAMs handle curved seasonality and sparse data", {
 set.seed(42);date<-seq(ts("2018-01-01"),by="month",length.out=60);t<-seq(0,5,length.out=60)
 dat<-tibble(calendar_month=date,n_icu_days=1000,n_detection_events=rnbinom(60,mu=20*exp(0.25*t+0.2*sin(2*pi*t)),size=10),n_tested=100,n_ns=rbinom(60,100,plogis(-2+0.4*t)))
 fit<-fit_temporal_model(dat,"n_icu_days");expect_equal(fit$summary$model_status,"estimated");expect_gt(fit$summary$annual_irr,1);expect_equal(nrow(fit$predictions),60L)
 fraction<-fit_temporal_model(dat,"n_tested","n_ns",TRUE);expect_equal(fraction$summary$model_status,"estimated");expect_true(all(fraction$predictions$fitted_value>=0 & fraction$predictions$fitted_value<=1));expect_gt(fraction$summary$annual_irr,1)
 expect_equal(fit_temporal_model(dat[1:5,],"n_icu_days")$summary$model_status,"insufficient_data")
})

test_that("CSV linkage identifiers preserve leading zeros", {
 f <- tempfile(fileext = ".csv"); readr::write_csv(tibble(patient_id = "001", hospitalization_id = "0002", organism_id = "00003"), f)
 x <- read_clif_csv(f, show_col_types = FALSE)
 expect_identical(x$organism_id, "00003"); expect_identical(x$hospitalization_id, "0002")
})

test_that("missing residual diagnostics cannot imply a confirmed increase", {
 d <- classify_temporal_direction(rep("estimated",3),rep(0.001,3),rep(20,3),c(NA,TRUE,FALSE))
 expect_equal(d,c("Residual diagnostic unavailable: exploratory","Residual dependence: exploratory","Increasing"))
})

test_that("bacterial influenzae names and generic bacilli are classified correctly", {
 expect_equal(classify_microbe_taxonomy(c("Haemophilus influenzae", "haemophilus_parainfluenzae", "influenza_a_virus", "gram_negative_bacillus")),c("Gram negative bacteria","Gram negative bacteria","Virus","Gram negative bacteria"))
})


test_that("season-adjusted curves match endpoint contrasts without changing inference", {
 set.seed(73)
 date <- seq(ts("2015-01-01"), by="month", length.out=120)
 t <- as.numeric(difftime(date,min(date),units="days"))/365.25
 seasonal <- sin(2*pi*(month(date)-1)/12)
 dat <- tibble(calendar_month=date,n_icu_days=1000,n_detection_events=rpois(120,100*exp(0.12*t+0.9*seasonal)))
 adjusted <- fit_temporal_model(dat,"n_icu_days")
 seasonal_fit <- fit_temporal_model(dat,"n_icu_days",season_adjusted=FALSE)
 expect_equal(adjusted$summary,seasonal_fit$summary)
 curve <- adjusted$predictions
 expect_equal((tail(curve$fitted_value,1)/head(curve$fitted_value,1))^(1/max(t)),adjusted$summary$annual_irr,tolerance=1e-7)
 expect_true(all(curve$ci_low <= curve$fitted_value & curve$ci_high >= curve$fitted_value))
 expect_lt(sd(diff(log(curve$fitted_value))),sd(diff(log(seasonal_fit$predictions$fitted_value)))/5)
})

test_that("mid-month carry-in never counts as a new admission", {
 stays<-tibble(icu_admission_month=ts("2020-01-01"),icu_in_dttm=ts("2020-01-05"),icu_in_dttm_clipped=ts("2020-01-15"),icu_out_dttm_clipped=ts("2020-01-25"))
 den<-monthly_icu_denominators(stays,ts("2020-01-01"))
 expect_equal(den$n_icu_admissions,0L);expect_equal(den$n_icu_days,10)
})

buddy_months<-seq(ts("2018-01-01"),by="month",length.out=48)
buddy_rows<-tibble(organism_id=paste0("o",seq_len(48)),culture_event_id=seq_len(48),patient_id="p",hospitalization_id="h",positive_culture=TRUE,organism_category=rep(c("escherichia_coli","staphylococcus_aureus"),each=24),fluid_category="blood_buffy",collect_dttm=buddy_months+days(1))
buddy_ast<-tibble(organism_id=paste0("o",seq_len(24)),antimicrobial_category="ceftriaxone",susceptibility_category="non_susceptible")
buddy_den<-tibble(calendar_month=buddy_months,n_icu_days=1000,n_icu_admissions=100L)
test_that("disappearance retains validated zero months but never invents tested fractions", {
 result<-build_susceptibility_analysis(buddy_rows,buddy_ast,buddy_den)
 x<-filter(result$monthly,specimen_stratum=="Overall")
 expect_equal(nrow(x),48L);expect_true(all(x$rate_model_eligible))
 expect_true(all(x$observation_status[25:48]=="organism_not_detected"))
 expect_true(all(x$non_susceptible_per_100_icu_days[25:48]==0))
 expect_true(all(is.na(x$non_susceptible_fraction[25:48])))
 gap<-build_susceptibility_analysis(buddy_rows[-48,],buddy_ast,buddy_den)$monthly
 expect_true(all(filter(gap,calendar_month==buddy_months[48])$observation_status=="source_unavailable"))
 incomplete<-buddy_rows;incomplete$organism_category[48]<-NA_character_
 unknown<-build_susceptibility_analysis(incomplete,buddy_ast,buddy_den)$monthly
 expect_true(all(filter(unknown,calendar_month==buddy_months[48])$observation_status=="organism_identification_incomplete"))
})
test_that("testing coverage includes isolates missing linkage IDs", {
 rows<-buddy_rows[rep(1,10),];rows$culture_event_id<-1:10;rows$organism_id<-c("a","b",rep(NA_character_,8))
 ast<-tibble(organism_id=c("a","b"),antimicrobial_category="ceftriaxone",susceptibility_category="susceptible")
 x<-build_susceptibility_analysis(rows,ast,buddy_den)$monthly %>% filter(calendar_month==buddy_months[1],specimen_stratum=="Overall")
 expect_equal(x$n_culture_isolates,10L);expect_equal(x$n_linkable_culture_isolates,2L)
 expect_equal(x$testing_fraction,0.2);expect_equal(x$testing_fraction_linkable,1)
 expect_equal(x$observation_status,"insufficient_linkage");expect_false(x$rate_model_eligible);expect_false(x$fraction_model_eligible)
})
test_that("AST availability distinguishes failed linkage and uninterpretable reports", {
 orphan<-buddy_ast;orphan$organism_id<-paste0("orphan",1:24)
 expect_equal(build_susceptibility_analysis(buddy_rows,orphan,buddy_den)$analysis_status,"no_linked_icu_tests")
 unknown<-buddy_ast;unknown$susceptibility_category<-NA_character_
 expect_equal(build_susceptibility_analysis(buddy_rows,unknown,buddy_den)$analysis_status,"no_interpretable_tests")
 empty<-filter(buddy_rows,FALSE)
 expect_equal(build_susceptibility_analysis(empty,buddy_ast,buddy_den)$analysis_status,"no_linkable_icu_isolates")
})

test_that("missing-ID fallback isolates with different source labels remain separate", {
 rows<-buddy_rows[rep(1,3),];rows$organism_id<-c("a",NA_character_,NA_character_);rows$organism_name<-c("E. coli","E. coli variant 1","E. coli variant 2")
 ast<-tibble(organism_id="a",antimicrobial_category="ceftriaxone",susceptibility_category="susceptible")
 x<-build_susceptibility_analysis(rows,ast,buddy_den)$monthly %>% filter(calendar_month==buddy_months[1],specimen_stratum=="Overall")
 expect_equal(x$n_culture_isolates,3L);expect_equal(x$testing_fraction,1/3)
})
