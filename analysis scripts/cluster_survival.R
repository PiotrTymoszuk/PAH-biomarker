## Checks for the differences in overall survival between the participant clusters

  insert_head()

# container list------
  
  cl_surv <- list()
  
# globals ----
  
  insert_msg('Globals setup')
  
  cl_surv$analysis_tbl <- clust[c('clust_obj_train', 
                                  'clust_obj_test')] %>% 
    map(extract, 'assignment') %>% 
    map(set_names, c('ID', 'clust_id')) %>% 
    map2(., multi_cox[c('score_IBK_0', 'score_LZ_0')], 
         left_join, by = 'ID') %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
  ## n numbers
  
  cl_surv$n_numbers <- cl_surv$analysis_tbl %>% 
    map(count, clust_id)

  cl_surv$n_tags <- cl_surv$n_numbers %>% 
    map(~map2_chr(.x$clust_id, .x$n, ~paste0(.x, ': n = ', .y))) %>% 
    map(paste, collapse = ', ') %>% 
    map(~paste0('\n', .x))
  
# Comparing the survival differences with Kaplan-Meier and Mentel-Henszel test -----
  
  insert_msg('Comparing the survial between the clusters')

  cl_surv$fits <- cl_surv$analysis_tbl %>% 
    map(~survfit(Surv(surv_months, death_study) ~ clust_id, 
                 data = .x))
  
  cl_surv$test_results <- cl_surv$analysis_tbl %>% 
    map(~survdiff(Surv(surv_months, death_study) ~ clust_id, 
                  data = .x, 
                  rho = 0))
  
  cl_surv$test_summary <- cl_surv$test_results %>% 
    map2_dfr(., names(.), ~tibble(cohort = .y, 
                                  chisq = .x$chisq, 
                                  df = length(.x$n) - 1)) %>% 
    mutate(p_value = 1 - pchisq(chisq, df), 
           p_adjusted = p.adjust(p_value, 'BH'))
  
# Plotting: Kaplan-Meier -----
  
  insert_msg('Drawing Kaplan-Meier plots')
  
  cl_surv$plots <- list(fit = cl_surv$fits, 
                        data = cl_surv$analysis_tbl, 
                        pval = signif(cl_surv$test_summary$p_adjusted, 2), 
                        title = globals$center_labs[c('IBK_0', 'LZ_0')]) %>% 
    pmap(ggsurvplot, 
         palette = unname(globals$cluster_colors[1:2]), 
         conf.int = TRUE, 
         xlab = 'Overall survival, months', 
         legend.title = 'Cluster', 
         legend.labs = c('#1', '#2'), 
         conf.int.alpha = 0.15, 
         pval.size = 2.75) %>% 
    map2(., cl_surv$n_tags, 
         ~.x$plot + 
           labs(tag = .y) + 
           globals$common_theme)
  
# END ----
  
  insert_tail()