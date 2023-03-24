## Checks for the differences in overall survival between the participant clusters

  insert_head()

# container list------
  
  cl_surv <- list()
  
# globals ----
  
  insert_msg('Globals setup')
  
  cl_surv$analysis_tbl <- clust$clust_obj %>% 
    map(extract, 'assignment') %>% 
    map(set_names, c('ID', 'clust_id')) %>% 
    map2(., multi_cox$lp_scores, 
         left_join, by = 'ID')
  
  ## n numbers

  cl_surv$clust_labs <- clust$clust_obj %>% 
    map(ngroups) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ': n = '))
  
  cl_surv$n_numbers <- pah_study[c("IBK_0", "LZ_0")] %>% 
    map(count, death_study)
  
  cl_surv$n_tags <- 
    list(paste0('total: n = ', sum(cl_surv$n_numbers$IBK_0$n), 
                ', events: n = ', cl_surv$n_numbers$IBK_0$n[2]), 
         paste0('total: n = ', sum(cl_surv$n_numbers$LZ_0$n), 
                ', events: n = ', cl_surv$n_numbers$LZ_0$n[2]))
  
# Comparing the survival differences with Kaplan-Meier and Mentel-Henszel test -----
  
  insert_msg('Comparing the survial between the clusters')

  cl_surv$fits <- cl_surv$analysis_tbl %>% 
    map(~surv_fit(Surv(surv_months, death_study) ~ clust_id, 
                  data = .x))

  cl_surv$test_summary <- cl_surv$fits %>% 
    map(surv_pvalue, method = 'survdiff') %>%
    compress(names_to = 'cohort') %>% 
    mutate(p_adjusted = p.adjust(pval, 'BH'), 
           significance = ifelse(p_adjusted >= 0.05, 
                                 paste0('ns (p = ', signif(p_adjusted, 2), ')'), 
                                 ifelse(p_adjusted < 0.001, 
                                        'p < 0.001', 
                                        paste('p =', signif(p_adjusted, 2)))))
  
# Plotting: Kaplan-Meier -----
  
  insert_msg('Drawing Kaplan-Meier plots')
  
  cl_surv$plots <- list(fit = cl_surv$fits, 
                        data = cl_surv$analysis_tbl, 
                        pval = cl_surv$test_summary$significance, 
                        title = globals$center_labs[c('IBK_0', 'LZ_0')], 
                        legend.labs = cl_surv$clust_labs) %>% 
    pmap(ggsurvplot, 
         palette = unname(globals$cluster_colors[1:2]), 
         conf.int = TRUE, 
         xlab = 'Overall survival, months', 
         legend.title = 'Cluster', 
         conf.int.alpha = 0.15, 
         pval.size = 2.75) %>% 
    map2(., cl_surv$n_tags, 
         ~.x$plot + 
           labs(subtitle = .y) + 
           globals$common_theme)
  
# END ----
  
  insert_tail()