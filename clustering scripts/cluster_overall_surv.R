# This script compares overall survival between the study clusters 

  insert_head()
  
# data container ----
  
  clust_os <- list()
  
# serial KM analysis -----
  
  insert_msg('KM analysis')
  
  ## KM objects
  
  clust_os$analysis_obj <- c('IBK_0', 'LZ_0') %>%
    map(function(cohort) model_km_cutpoint(inp_table = study_clust$analysis_tbl[[cohort]], 
                                           surv_object = pah_study$surv_obj[[cohort]], 
                                           indep_variable = 'clust_id')) %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
  ## KM summary
  
  clust_os$summaries <- clust_os$analysis_obj %>% 
    map_dfr(km_summary) %>% 
    mutate(cohort = names(clust_os$analysis_obj))
  
  ## plots
  
  clust_os$plots <- list(km_model_object = clust_os$analysis_obj, 
                         plot_title = globals$center_labs[c('IBK_0', 'LZ_0')]) %>% 
    pmap(plot_km, 
         x_lab = 'OS, months', 
         palette = unname(study_clust$clust_colors), 
         legend.labs = study_clust$clust_labs)
# END ----
  
  insert_tail()


  
