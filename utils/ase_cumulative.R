# Observed cumulative detection/collection, with fixed new-admission denominators.
# Upstream culture-core filtering enforces study-window, half-open ICU intervals.
build_ase_cumulative <- function(z, max_day=14L, max_hour=168L) {
  if(any(!is.finite(c(max_day,max_hour))) || any(c(max_day,max_hour)<1) || any(c(max_day,max_hour)!=floor(c(max_day,max_hour)))) stop("Cumulative timing horizons must be positive integers")
  stays <- z$stays %>% filter(icu_in_dttm>=icu_in_dttm_clipped)
  if(anyDuplicated(stays$icu_admission_id))stop("Cumulative curves require unique ICU admissions")
  den <- tibble(ase_group=ase_group_levels) %>% left_join(count(stays,ase_group,name="n_icu_admissions"),by="ase_group") %>% mutate(n_icu_admissions=coalesce(n_icu_admissions,0L))
  observation <- function(hours) bind_rows(lapply(ase_group_levels,function(g) {
    durations <- filter(stays,ase_group==g)$icu_los_days*24
    tibble(ase_group=g,icu_hour=hours,n_icu_admissions_under_observation=vapply(hours,function(h)sum(durations>h),integer(1)))
  })) %>% left_join(den,by="ase_group")
  obs <- observation(0:max_hour)
  day_obs <- observation(24*(0:max_day)) %>% mutate(icu_day=icu_hour/24)
  events <- z$events %>% semi_join(stays,by="icu_admission_id") %>%
    mutate(hours=as.numeric(difftime(collect_dttm,icu_in_dttm,units="hours")),event_hour=ceiling(hours)) %>% filter(is.finite(hours),hours>=0)
  # At each day endpoint count a first event at or before the exact endpoint.
  first <- bind_rows(events %>% group_by(ase_group,icu_admission_id) %>% summarise(hours=if(n())min(hours) else NA_real_,.groups="drop") %>% mutate(outcome="First culture"),
    events %>% filter(any_positive_culture) %>% group_by(ase_group,icu_admission_id) %>% summarise(hours=if(n())min(hours) else NA_real_,.groups="drop") %>% mutate(outcome="First positive culture")) %>%
    mutate(icu_day=ceiling(hours/24)) %>% filter(icu_day<=max_day) %>% count(ase_group,outcome,icu_day,name="n_first_events")
  first <- tidyr::crossing(ase_group=ase_group_levels,outcome=c("First culture","First positive culture"),icu_day=0:max_day) %>%
    left_join(first,by=c("ase_group","outcome","icu_day")) %>% mutate(n_first_events=coalesce(n_first_events,0L)) %>%
    group_by(ase_group,outcome) %>% arrange(icu_day,.by_group=TRUE) %>% mutate(n_admissions_with_event_by_day=cumsum(n_first_events)) %>% ungroup() %>%
    left_join(select(day_obs,-icu_hour),by=c("ase_group","icu_day")) %>%
    mutate(cumulative_percent=if_else(n_icu_admissions>0,100*n_admissions_with_event_by_day/n_icu_admissions,NA_real_))
  rows <- z$rows %>% semi_join(stays,by="icu_admission_id") %>% filter(positive_culture,!is.na(organism_category),!organism_category %in% c("na","unknown","missing")) %>%
    mutate(hours=as.numeric(difftime(collect_dttm,icu_in_dttm,units="hours"))) %>% filter(is.finite(hours),hours>=0)
  totals <- rows %>% distinct(ase_group,culture_event_id,organism_category) %>% count(ase_group,organism_category,name="n_total_detection_events")
  org <- rows %>% group_by(ase_group,icu_admission_id,organism_category) %>% summarise(icu_hour=if(n())ceiling(min(hours)) else NA_real_,.groups="drop") %>%
    filter(icu_hour<=max_hour) %>% count(ase_group,organism_category,icu_hour,name="n_first_detections")
  org <- tidyr::crossing(ase_group=ase_group_levels,organism_category=sort(unique(rows$organism_category)),icu_hour=0:max_hour) %>%
    left_join(org,by=c("ase_group","organism_category","icu_hour")) %>% mutate(n_first_detections=coalesce(n_first_detections,0L)) %>%
    group_by(ase_group,organism_category) %>% arrange(icu_hour,.by_group=TRUE) %>% mutate(n_admissions_with_organism_by_hour=cumsum(n_first_detections)) %>% ungroup() %>%
    left_join(obs,by=c("ase_group","icu_hour")) %>% left_join(totals,by=c("ase_group","organism_category")) %>%
    mutate(n_total_detection_events=coalesce(n_total_detection_events,0L),cumulative_percent=if_else(n_icu_admissions>0,100*n_admissions_with_organism_by_hour/n_icu_admissions,NA_real_))
  totals <- events %>% mutate(fluid_category=coalesce(fluid_category,"missing")) %>% count(ase_group,fluid_category,name="n_total_culture_events")
  cultures <- events %>% mutate(fluid_category=coalesce(fluid_category,"missing")) %>% filter(event_hour<=max_hour) %>% count(ase_group,fluid_category,icu_hour=event_hour,name="n_culture_events")
  cultures <- tidyr::crossing(ase_group=ase_group_levels,fluid_category=sort(unique(totals$fluid_category)),icu_hour=0:max_hour) %>%
    left_join(cultures,by=c("ase_group","fluid_category","icu_hour")) %>% mutate(n_culture_events=coalesce(n_culture_events,0L)) %>%
    group_by(ase_group,fluid_category) %>% arrange(icu_hour,.by_group=TRUE) %>% mutate(n_culture_events_by_hour=cumsum(n_culture_events)) %>% ungroup() %>%
    left_join(obs,by=c("ase_group","icu_hour")) %>% left_join(totals,by=c("ase_group","fluid_category")) %>%
    mutate(n_total_culture_events=coalesce(n_total_culture_events,0L),cumulative_events_per_100_admissions=if_else(n_icu_admissions>0,100*n_culture_events_by_hour/n_icu_admissions,NA_real_))
  qc <- den %>% mutate(n_carry_in_stays_excluded=vapply(ase_group,function(g)sum(z$stays$ase_group==g & z$stays$icu_in_dttm<z$stays$icu_in_dttm_clipped),integer(1)),
    max_icu_day=max_day,max_icu_hour=max_hour,denominator="all subgroup ICU admissions starting in study window",estimator="observed cumulative proportion; no censoring or competing-risk adjustment")
  list(first=first,organisms=org,cultures=cultures,observation=obs,qc=qc)
}
