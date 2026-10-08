# Timing distributions from privacy-safe hourly culture counts; includes repeat events.
build_ase_culture_density <- function(cultures, bandwidth=6, bin_hours=6) {
  horizon<-max(cultures$icu_hour)
  if(!is.finite(horizon)||horizon<=0||!is.finite(bandwidth)||bandwidth<=0)stop("Invalid culture density horizon or bandwidth")
  hourly<-bind_rows(
    cultures %>% group_by(ase_group,icu_hour) %>% summarise(n_culture_events=sum(n_culture_events),.groups="drop") %>% mutate(specimen="All cultures"),
    cultures %>% filter(fluid_category=="blood_buffy") %>% select(ase_group,icu_hour,n_culture_events) %>% mutate(specimen="Blood cultures"))
  grid<-seq(0,horizon,length.out=673)
  curves<-list();bins<-list();hourly_bins<-list();k<-0L
  for(spec in c("All cultures","Blood cultures"))for(g in ase_group_levels) {
    x<-hourly %>% filter(specimen==spec,ase_group==g)
    n<-sum(x$n_culture_events)
    # Upstream bins use ceiling(time). Approximate nonzero bin times by their midpoints.
    times<-pmax(0,x$icu_hour-.5)
    y<-rep(NA_real_,length(grid))
    if(n>0) {
      y<-vapply(grid,function(t)sum(x$n_culture_events*(dnorm((t-times)/bandwidth)+dnorm((t+times)/bandwidth)+dnorm((t-(2*horizon-times))/bandwidth)))/(n*bandwidth),numeric(1))
      area<-sum(diff(grid)*(head(y,-1)+tail(y,-1))/2)
      y<-y/area
    }
    k<-k+1L
    curves[[k]]<-tibble(ase_group=g,specimen=spec,icu_hour=grid,density_per_hour=y,n_culture_events_within_horizon=n,bandwidth_hours=bandwidth,max_icu_hour=horizon)
    b<-tibble(bin_end_hour=pmin(seq(bin_hours,horizon+bin_hours,by=bin_hours),horizon)) %>% distinct() %>% mutate(bin_start_hour=lag(bin_end_hour,default=0))
    counts<-x %>% mutate(bin_end_hour=pmin(pmax(bin_hours,ceiling(icu_hour/bin_hours)*bin_hours),horizon)) %>% group_by(bin_end_hour) %>% summarise(n_culture_events=sum(n_culture_events),.groups="drop")
    bins[[k]]<-b %>% left_join(counts,by="bin_end_hour") %>% mutate(n_culture_events=coalesce(n_culture_events,0L),ase_group=g,specimen=spec,n_culture_events_within_horizon=n)
    hourly_counts<-x %>% mutate(bin_end_hour=pmax(1,icu_hour)) %>% group_by(bin_end_hour) %>% summarise(n_culture_events=sum(n_culture_events),.groups="drop")
    hourly_bins[[k]]<-tibble(bin_end_hour=seq_len(horizon),bin_start_hour=seq_len(horizon)-1) %>% left_join(hourly_counts,by="bin_end_hour") %>% mutate(n_culture_events=coalesce(n_culture_events,0L),ase_group=g,specimen=spec,n_culture_events_within_horizon=n)
  }
  specimen_counts<-cultures %>% mutate(bin_end_hour=pmax(1,icu_hour)) %>% group_by(ase_group,fluid_category,bin_end_hour) %>% summarise(n_culture_events=sum(n_culture_events),.groups="drop") %>% mutate(bin_start_hour=bin_end_hour-1)
  list(curves=bind_rows(curves),counts=bind_rows(bins),hourly_counts=bind_rows(hourly_bins),specimen_hourly_counts=specimen_counts)
}

plot_ase_culture_density <- function(x, out, site, stamp) {
  colors<-setNames(c("#0072B2","#E69F00","#CC79A7"),ase_group_levels)
  labels<-function(x)stringr::str_wrap(x,24)
  horizon<-max(x$curves$max_icu_hour)
  caption<-paste0("All collections, including repeats; new ICU admissions only. Groups reflect the full hospitalization.\nDensity integrates to 1 within 0–",horizon," hours for each group/specimen; common ",first(x$curves$bandwidth_hours),"-hour bandwidth, reflected boundaries. Hourly counts approximate collection times.")
  p<-ggplot(filter(x$curves,n_culture_events_within_horizon>0),aes(icu_hour,density_per_hour,color=ase_group,fill=ase_group))+geom_area(alpha=.12,position="identity")+geom_line(linewidth=.8)+facet_wrap(~specimen,ncol=1)+scale_color_manual(values=colors,limits=ase_group_levels,labels=labels,name=NULL)+scale_fill_manual(values=colors,limits=ase_group_levels,labels=labels,name=NULL)+labs(x="Hours since ICU admission",y="Culture timing density (per hour)",caption=NULL)+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,paste0("ase_culture_timing_density_",site,"_",stamp,".png")),p,width=12,height=8,dpi=200)
  p<-ggplot(x$counts,aes((bin_start_hour+bin_end_hour)/2,n_culture_events,color=ase_group))+geom_line(linewidth=.8)+geom_point(size=1)+facet_wrap(~specimen,ncol=1,scales="free_y")+scale_color_manual(values=colors,limits=ase_group_levels,labels=labels,name=NULL)+labs(x="Hours since ICU admission (6-hour bin midpoints)",y="Number of culture events",caption=NULL)+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,paste0("ase_culture_timing_counts_",site,"_",stamp,".png")),p,width=12,height=8,dpi=200)
  p<-ggplot(x$hourly_counts,aes(xmin=bin_start_hour,xmax=bin_end_hour,ymin=0,ymax=n_culture_events,color=ase_group,fill=ase_group))+geom_rect(alpha=.22,linewidth=.25,position="identity")+facet_wrap(~specimen,ncol=1,scales="free_y")+scale_color_manual(values=colors,limits=ase_group_levels,labels=labels,name=NULL)+scale_fill_manual(values=colors,limits=ase_group_levels,labels=labels,name=NULL)+scale_x_continuous(breaks=seq(0,horizon,by=24),expand=expansion(mult=c(0,.01)))+labs(x="Hours since ICU admission",y="Number of culture events per 1-hour bin",caption=NULL)+theme_bw()+theme(legend.position="bottom")
  ggsave(file.path(out,paste0("ase_culture_timing_histogram_1h_",site,"_",stamp,".png")),p,width=12,height=8,dpi=200)
}

# Every nonempty non-blood standardized source; retain zero-volume sources in numeric exports.
plot_ase_specimen_histograms <- function(x, out, site, stamp) {
  colors<-setNames(c("#0072B2","#E69F00","#CC79A7"),ase_group_levels)
  sources<-x$specimen_hourly_counts %>% group_by(fluid_category) %>% summarise(n=sum(n_culture_events),.groups="drop") %>% filter(n>0,fluid_category!="blood_buffy") %>% arrange(desc(n),fluid_category) %>% pull(fluid_category)
  if(!length(sources))return(invisible(NULL))
  for(page in seq_len(ceiling(length(sources)/9))) {
    selected<-sources[seq((page-1)*9+1,min(page*9,length(sources)))]
    dat<-x$specimen_hourly_counts %>% filter(fluid_category %in% selected) %>% mutate(fluid_category=factor(fluid_category,levels=selected),ase_group=factor(ase_group,levels=ase_group_levels))
    p<-ggplot(dat,aes(xmin=bin_start_hour,xmax=bin_end_hour,ymin=0,ymax=n_culture_events,color=ase_group,fill=ase_group))+geom_rect(alpha=.22,linewidth=.2,position="identity")+facet_wrap(~fluid_category,ncol=3,scales="free_y",labeller=labeller(fluid_category=function(v)stringr::str_wrap(gsub("_"," ",v),30)))+scale_color_manual(values=colors,limits=ase_group_levels,labels=function(v)stringr::str_wrap(v,24),name=NULL)+scale_fill_manual(values=colors,limits=ase_group_levels,labels=function(v)stringr::str_wrap(v,24),name=NULL)+scale_x_continuous(breaks=c(0,48,96,144),expand=expansion(mult=c(0,.01)))+labs(x="Hours since ICU admission",y="Culture events per 1-hour bin",caption=NULL)+theme_bw()+theme(legend.position="bottom")
    ggsave(file.path(out,paste0("ase_other_specimen_histogram_1h_",site,"_",stamp,"_page",page,".png")),p,width=14,height=10,dpi=200)
  }
}
