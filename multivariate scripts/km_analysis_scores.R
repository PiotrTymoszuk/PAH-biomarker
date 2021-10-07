# This script analyses differences in survival in the study cohorts in respect to the developed
# signature score values. The scores values are stratified by quartiles

  insert_head()
  
# data container -----
  
  km <- list()
  
# globals: variables and analysi tables -----
  
  ## variables
  
  km$sign_scores <- multi_modeling$cv_pass$stats$model_id
  km$comparators <- pah_study$comparators$variable
  km$variables <- c(km$sign_scores, 
                    km$comparators)
  
  ## colors
  
  km$cmm_palette <- c('steelblue', 
                      'cornsilk4', 
                      'coral3', 
                      'firebrick4', 
                      'black')
  
  ## analysis table with the sgnature scores stratified by quartiles separately for each cohort
  
  km$analysis_tbl <- multi_modeling$score_tbl %>% 
    map(function(cohort) cohort %>% 
          map_dfc(cut_quartile)) %>% 
    map2(., 
         map(pah_study[c('IBK_0', 'LZ_0')], 
             select, 
             all_of(km$comparators)), 
         cbind)
  
# serial KM analysis ----
  
  insert_msg('Serial KM analysis')
  
  ## km objects
  
  km$analysis_obj <- c('IBK_0', 'LZ_0') %>% 
    map(function(cohort) km$variables %>% 
          map(safely(model_km_cutpoint), 
              inp_table = km$analysis_tbl[[cohort]], 
              surv_object = pah_study$surv_obj[[cohort]])) %>% 
    map(function(cohort) cohort %>% 
          map(~.x$result) %>% 
          set_names(km$variables) %>% 
          compact) %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
  ## km summaries
  
  km$summaries <- km$analysis_obj %>% 
    map(km_summary) %>% 
    map(mutate, 
        p_adj = p.adjust(P_value, 'BH'))
  
  ## KM plots
  
  km$plots <- c('IBK_0', 'LZ_0') %>% 
    map(function(cohort) list(km_model_object = km$analysis_obj[[cohort]], 
                              p_value = signif(km$summaries[[cohort]][['p_adj']], 2), 
                              plot_title = ifelse(names(km$analysis_obj[[cohort]]) %in% pah_study$comparators$variable, 
                                                  globals$comp_labs[names(km$analysis_obj[[cohort]])], 
                                                  stri_replace(names(km$analysis_obj[[cohort]]), 
                                                               fixed = 'sign_', 
                                                               replacement = 'Signature ')) %>% 
                                                  paste(., globals$center_labs[[cohort]], sep = ': '), 
                              legend.labs = map(km$analysis_obj[[cohort]], 
                                                ~.x$test_table$var_strata) %>% 
                                map(factor) %>% 
                                map(levels)) %>% 
          pmap(plot_km, 
               palette = km$cmm_palette, 
               x_lab = 'OS, months', 
               legend.title = '') %>% 
          set_names(names(km$analysis_obj[[cohort]]))) %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
# END ----
  
  insert_tail()