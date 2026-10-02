suppressPackageStartupMessages({library(dplyr);library(tidyr);library(lubridate);library(testthat);library(readr)})
source("utils/susceptibility.R");source("utils/trends.R");source("utils/pooling.R")
set.seed(160)
months <- seq(as.Date("2018-01-01"),by="month",length.out=60)
fixture <- crossing(site_name=paste0("S",1:4),calendar_month=months) %>% mutate(site_index=as.integer(sub("S","",site_name)),t=as.numeric(calendar_month-min(calendar_month))/365.25,n_icu_days=1000*site_index,n_icu_admissions=100*site_index,n_observed_culture_events=100L,n_detection_events=rnbinom(n(),mu=n_icu_days*0.06*exp(0.15*t+0.25*site_index+0.04*(site_index-2.5)*(t-2)^2+0.1*site_index*sin(2*pi*(month(calendar_month)-1)/12)),size=100),organism_category="escherichia_coli")
test_that("common months and validation prevent incompatible contrasts", {
 dat <- fixture %>% filter(!(site_name=="S2" & calendar_month==months[10]),!(site_name=="S3" & calendar_month<months[13]))
 x <- common_pooling_months(dat,"n_icu_days","n_detection_events")
 expect_equal(n_distinct(x$calendar_month),48L);expect_true(all(table(x$site_name)==48))
 dat <- fixture;dat$n_observed_culture_events[1]<-0
 expect_equal(n_distinct(common_pooling_months(dat,"n_icu_days","n_detection_events")$calendar_month),59L)
 expect_error(common_pooling_months(bind_rows(fixture,fixture[1,]),"n_icu_days","n_detection_events"),"Duplicate")
 expect_error(pooling_forbidden(tibble(patient_id="x")),"private")
})
test_that("site models are refit over common months and pool without significance selection", {
 result <- pool_screen(fixture,"n_icu_days","n_detection_events")
 expect_true(all(result$effects$model_status=="estimated"));expect_equal(result$meta$n_sites,4L)
 expect_equal(result$meta$model_status,"estimated");expect_gt(result$meta$annual_effect,1)
 expect_true(all(is.finite(result$effects$theta)));expect_true(all(is.finite(result$effects$long_term_edf)))
 effects <- result$effects;effects$residual_dependence_flag[1]<-TRUE
 expect_equal(meta_pool(effects)$inference_status,"exploratory_residual_dependence")
 expect_equal(meta_pool(effects[1,])$model_status,"insufficient_sites")
 effects$log_annual_effect<-c(-0.001,0.001,0.002,-0.002);effects$log_annual_se<-0.1
 expect_equal(meta_pool(effects)$n_sites,4L)
})
test_that("joint curves average sites equally and include uncertainty and diagnostics", {
 result <- pool_screen(fixture,"n_icu_days","n_detection_events")
 expect_equal(result$joint$model_status,"estimated")
 curves<-result$curves;means<-curves %>% filter(curve_type=="site") %>% group_by(calendar_month) %>% summarise(mean_value=mean(fitted_value),.groups="drop")
 pooled<-curves %>% filter(curve_type=="pooled_equal_site") %>% left_join(means,by="calendar_month")
 expect_equal(pooled$fitted_value,pooled$mean_value,tolerance=1e-10)
 expect_true(all(pooled$ci_low<=pooled$fitted_value & pooled$ci_high>=pooled$fitted_value))
 expect_equal(nrow(result$diagnostics),4L)
 weighted<-curves %>% filter(curve_type=="site") %>% mutate(weight=as.numeric(sub("S","",site_name))) %>% group_by(calendar_month) %>% summarise(weighted_mean=weighted.mean(fitted_value,weight),.groups="drop")
 expect_gt(max(abs(pooled$fitted_value-weighted$weighted_mean)),0.1)
 expect_equal(pool_screen(filter(fixture,site_name %in% c("S1","S2")),"n_icu_days","n_detection_events")$joint$model_status,"insufficient_sites")
})
test_that("AST uses shared testing coverage and tested denominators", {
 ast <- fixture %>% mutate(n_interpretable=200L,n_non_susceptible=rbinom(n(),200,plogis(-2+0.2*t)),n_susceptible=n_interpretable-n_non_susceptible,testing_fraction=0.9,n_culture_isolates=222L,n_linkable_culture_isolates=222L,n_positive_rows_missing_organism_category=0L,culture_coverage_validated=TRUE)
 ast$n_culture_isolates[1]<-ast$n_linkable_culture_isolates[1]<-2000L
 rate<-pool_screen(ast,"n_icu_days","n_non_susceptible")
 expect_true(all(rate$effects$n_common_months==59))
 fraction<-pool_screen(ast,"n_interpretable","n_non_susceptible",TRUE)
 expect_true(all(fraction$effects$n_common_months==60))
 expect_equal(fraction$joint$model_status,"estimated")
 expect_true(all(fraction$curves$fitted_value>0 & fraction$curves$fitted_value<1))
 expect_error(common_pooling_months(mutate(ast,n_non_susceptible=201),"n_interpretable","n_non_susceptible",TRUE),"exceeds")
})

test_that("pooled AST preserves validated zeros and excludes poor linkage", {
 ast <- fixture %>% mutate(n_interpretable=20L,n_non_susceptible=10L,n_culture_isolates=20L,n_linkable_culture_isolates=20L,n_positive_rows_missing_organism_category=0L,culture_coverage_validated=TRUE,testing_fraction=1)
 ast[1,c("n_interpretable","n_non_susceptible","n_culture_isolates","n_linkable_culture_isolates")] <- 0L
 expect_equal(n_distinct(common_pooling_months(ast,"n_icu_days","n_non_susceptible")$calendar_month),60L)
 expect_equal(n_distinct(common_pooling_months(ast,"n_interpretable","n_non_susceptible",TRUE)$calendar_month),59L)
 ast$culture_coverage_validated[1]<-FALSE
 expect_equal(n_distinct(common_pooling_months(ast,"n_icu_days","n_non_susceptible")$calendar_month),59L)
 ast$n_culture_isolates[1]<-100L;ast$n_linkable_culture_isolates[1]<-20L;ast$n_interpretable[1]<-20L
 expect_equal(n_distinct(common_pooling_months(ast,"n_interpretable","n_non_susceptible",TRUE)$calendar_month),59L)
})

# Aggregate fixture runs exercise registry, provenance, CLI, figures and optional AST.
scratch <- Sys.getenv("POOLING_TEST_OUTPUT",tempfile("clif_pooling_"));
if(dir.exists(scratch))stop("Use a fresh POOLING_TEST_OUTPUT directory.")
dir.create(scratch,recursive=TRUE)
registry <- tibble(site_name=paste0("S",1:4),run_dir=file.path(scratch,paste0("S",1:4)),validated_start_date=as.Date("2018-01-01"),validated_end_date=as.Date("2022-12-31"),culture_qc_pass=TRUE,ast_qc_pass=c(TRUE,TRUE,TRUE,FALSE))
source_files<-c("utils/culture_core.R","utils/trends.R","utils/susceptibility.R","code/08_organism_trends.R","code/09_susceptibility_trends.R","config/mcide/clif_microbiology_susceptibility_category.csv","config/mcide/clif_microbiology_susceptibility_antibiotics_category.csv")
for(i in 1:4){
 r<-registry[i,];for(folder in c("provenance","site_exports","organism_trends","susceptibility"))dir.create(file.path(r$run_dir,folder),recursive=TRUE)
 jsonlite::write_json(list(site_name=r$site_name,run_id=paste0("synthetic_",i),analysis_status="completed",study_start_date="2018-01-01",study_end_date="2022-12-31",code_and_mapping_md5=as.list(tools::md5sum(source_files))),file.path(r$run_dir,"provenance","run_manifest.json"),auto_unbox=TRUE)
 dat<-filter(fixture,site_name==r$site_name)
 write_csv(select(dat,calendar_month,n_icu_days,n_icu_admissions,n_observed_culture_events),file.path(r$run_dir,"organism_trends","monthly_icu_denominators_for_organism_trends_fixture.csv"))
 write_csv(select(dat,calendar_month,organism_category,n_detection_events,site_name),file.path(r$run_dir,"organism_trends","monthly_all_organism_detection_counts_fixture.csv"))
 write_csv(tibble(file=character(),column=character()),file.path(r$run_dir,"site_exports","site_export_privacy_audit_fixture.csv"))
 if(i<4){
   ast<-dat %>% mutate(antimicrobial_category="ceftriaxone",specimen_stratum="Overall",n_interpretable=200L,n_non_susceptible=rbinom(n(),200,plogis(-2+0.2*t)),n_susceptible=n_interpretable-n_non_susceptible,testing_fraction=0.9,n_culture_isolates=222L,n_linkable_culture_isolates=222L,n_positive_rows_missing_organism_category=0L,culture_coverage_validated=TRUE) %>% select(site_name,calendar_month,organism_category,antimicrobial_category,specimen_stratum,n_interpretable,n_non_susceptible,n_susceptible,testing_fraction,n_icu_days,n_culture_isolates,n_linkable_culture_isolates,n_observed_culture_events,n_positive_rows_missing_organism_category,culture_coverage_validated)
   write_csv(ast,file.path(r$run_dir,"susceptibility","monthly_organism_antimicrobial_susceptibility_fixture.csv"))
 }
}
test_that("registry checks complete runs, valid intervals and all-category zero completion", {
 inputs<-read_pooling_sites(registry);expect_equal(nrow(inputs$culture),240L);expect_equal(n_distinct(inputs$ast$site_name),3L)
 expect_error(read_pooling_sites(mutate(registry,validated_start_date=as.Date("2017-01-01"))),"outside")
 expect_error(read_pooling_sites(bind_rows(registry,registry[1,])),"duplicate")
 # A category absent from one ALL-organism export is correctly completed as zero.
 f<-file.path(registry$run_dir[4],"organism_trends","monthly_all_organism_detection_counts_fixture.csv")
 original<-read_csv(f,show_col_types=FALSE);write_csv(original %>% mutate(organism_category="candida_albicans"),f)
 inputs<-read_pooling_sites(registry)
 expect_true(all(filter(inputs$culture,site_name=="S4",organism_category=="escherichia_coli")$n_detection_events==0))
 write_csv(original,f)
 # Partial calendar months are dropped.
 partial<-registry;partial$validated_start_date[1]<-as.Date("2018-01-15")
 expect_equal(min(filter(read_pooling_sites(partial)$culture,site_name=="S1")$calendar_month),as.Date("2018-02-01"))
})
test_that("registry rejects wrong-site AST, failed provenance and incompatible source", {
 f<-file.path(registry$run_dir[1],"susceptibility","monthly_organism_antimicrobial_susceptibility_fixture.csv")
 original<-read_csv(f,show_col_types=FALSE);write_csv(mutate(original,site_name="WRONG_SITE"),f)
 expect_error(read_pooling_sites(registry),"AST export site")
 write_csv(original,f)
 manifest_path<-file.path(registry$run_dir[1],"provenance","run_manifest.json")
 original_manifest<-jsonlite::fromJSON(manifest_path,simplifyVector=FALSE)
 failed<-original_manifest;failed$analysis_status<-"failed"
 jsonlite::write_json(failed,manifest_path,auto_unbox=TRUE)
 expect_error(read_pooling_sites(registry),"incomplete")
 incompatible<-original_manifest;incompatible$code_and_mapping_md5[["utils/trends.R"]]<-"changed"
 jsonlite::write_json(incompatible,manifest_path,auto_unbox=TRUE)
 expect_error(read_pooling_sites(registry),"Central and site")
 jsonlite::write_json(original_manifest,manifest_path,auto_unbox=TRUE)
})
registry_path<-file.path(scratch,"registry.csv");write_csv(registry,registry_path)
out<-file.path(scratch,"pooled");log<-file.path(scratch,"pooling.log")
status<-system2(file.path(R.home("bin"),"Rscript"),c("code/11_pool_site_trends.R",shQuote(registry_path),shQuote(out)),stdout=log,stderr=log)
if(status!=0){cat(readLines(log),sep="\n");stop("Synthetic central pipeline failed.")}
test_that("central CLI writes complete pooled outputs and figures", {
 expect_equal(jsonlite::fromJSON(file.path(out,"pooling_manifest.json"))$analysis_status,"completed")
 meta<-read_csv(file.path(out,"pooled_meta.csv"),show_col_types=FALSE)
 expect_equal(nrow(meta),5L);expect_true(all(meta$model_status=="estimated"));expect_true(all(is.finite(meta$fdr_p_value)))
 expect_equal(length(list.files(out,pattern="_forest[.]png$")),5L);expect_equal(length(list.files(out,pattern="_joint[.]png$")),5L)
 expect_true(file.exists(file.path(out,"input_file_hashes.csv")))
})
cat("Synthetic four-site pooling passed. Artifacts: ",out,"\n",sep="")
