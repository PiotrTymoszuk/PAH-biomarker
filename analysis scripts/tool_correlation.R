# Correlation of the risk tool values in the study cohorts
#
# The purpose is to check to which extent the established risk assessment tools
# overlap with each other and with the newly established Elastic net score
#
# The Reveal scores are skipped from the final analysis

  insert_head()
  
# container -----
  
  corr_tools <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## analysis variables
  
  corr_tools$var_lexicon <- 
    rbind(tibble(variable = 'score', 
                 label = 'ElasticNet'), 
          pah_study$comparators[c('variable', 'label')]) %>% 
    filter(!stri_detect(variable, fixed = 'Reveal'))
  
  ## analysis tables
  
  corr_tools$analysis_tbl <- 
    map2(multi_cox$lp_scores, 
         pah_study$data_master %>% 
           blast(timepoint) %>% 
           map(select, 
               ID, any_of(corr_tools$var_lexicon$variable)), 
         left_join, by = 'ID') %>% 
    map(~map_dfc(.x, function(x) if(is.factor(x)) as.numeric(x) else x)) %>% 
    map(column_to_rownames, 'ID')
  
  #corr_tools$analysis_tbl$LZ_0 <-  corr_tools$analysis_tbl$LZ_0 %>% 
   # select(-starts_with('Reveal'))
  
  ## cohort specific variables (no REVEAL for the LZ/W cohort!)
  
  corr_tools$variables <- corr_tools$analysis_tbl %>% 
    map(select, any_of(corr_tools$var_lexicon$variable)) %>% 
    map(names)
  
  corr_tools$analysis_tbl <- 
    map2(corr_tools$analysis_tbl, 
         corr_tools$variables, 
         ~select(.x, all_of(.y)))
  
  ## median centering of the tables
  
  corr_tools$analysis_tbl <- corr_tools$analysis_tbl %>% 
    map(center_data, 'median')
  
  ## variable pairs
  
  corr_tools$var_pairs <- corr_tools$variables %>% 
    map(~combn(x = .x, m = 2, simplify = FALSE))
  
# Pairwise correlations -------
  
  insert_msg('Pairwise correlations')
  
  ## Spearman's rank test
  
  corr_tools$test <- 
    map2(corr_tools$analysis_tbl,
         corr_tools$var_pairs, 
         function(data, pair) pair %>% 
           map_dfr(~correlate_variables(data, 
                                        variables = .x, 
                                        what = 'correlation', 
                                        type = 'spearman', 
                                        pub_styled = FALSE))) %>% 
    map(mutate, 
        p_adjusted = p.adjust(p_value, 'BH'), 
        significance = ifelse(p_adjusted >= 0.05, 
                              paste0('ns (p = ', signif(p_adjusted, 2), ')'), 
                              ifelse(p_adjusted < 0.001, 
                                     'p < 0.001', 
                                     paste('p =', signif(p_adjusted, 2)))))
  
# Kendall's concordance coefficient ------
  
  insert_msg('KCC')
  
  corr_tools$kcc <- corr_tools$analysis_tbl %>% 
    map(as.matrix) %>% 
    map(KendallW, 
        correct = TRUE, 
        test = TRUE) %>% 
    map(~.x[c('estimate', 'statistic', 'p.value')]) %>% 
    map(as_tibble) %>% 
    compress(names_to = 'cohort')
  
# Correlograms --------
  
  insert_msg('Correlograms')
  
  corr_tools$plots <- 
    list(x = corr_tools$test, 
         v = corr_tools$variables, 
         y = c('Training: IBK', 'Test: LZ/W'), 
         z = paste('KCC = ', signif(corr_tools$kcc$estimate, 2))) %>% 
    pmap(function(x, v, y, z) x %>% 
           ggplot(aes(x = variable1,
                      y = variable2, 
                      fill = estimate, 
                      size = abs(estimate))) + 
           geom_point(shape = 21) + 
           geom_text(aes(label = signif(estimate, 2)), 
                     size = 2.5, 
                     hjust = 0.5,
                     vjust = -1.4) + 
           scale_x_discrete(limits = v, 
                            labels = exchange(v, 
                                              dict = corr_tools$var_lexicon)) + 
           scale_y_discrete(limits = v, 
                            labels = exchange(v, 
                                              dict = corr_tools$var_lexicon)) + 
           scale_fill_gradient2(low = 'steelblue',
                                mid = 'white', 
                                high = 'firebrick', 
                                midpoint = 0.5, 
                                limits = c(0, 1), 
                                name = expression(rho)) + 
           scale_radius(limits = c(0, 1), 
                        range = c(0.5, 4.5), 
                        name = expression(rho)) + 
           guides(fill = 'legend', 
                  size = 'legend', 
                  alpha = 'none') + 
           globals$common_theme + 
           theme(axis.title = element_blank()) + 
           labs(title = y, 
                subtitle = z))
  
# Factor analysis ------
  
  insert_msg('Factor analysis')
  
  ## ignoring the reveal-family variables absent from the test cohort
  
  ## with three components 
  ## see the output of omega
  
  corr_tools$fa_obj <- corr_tools$analysis_tbl %>% 
    map(select, -starts_with('Reveal')) %>% 
    map(reduce_data, 
        kdim = 3, 
        red_fun = 'fa')
  
  ## McDonald's omega
  
  corr_tools$omega <- corr_tools$analysis_tbl %>% 
    map(select, -starts_with('Reveal')) %>% 
    map(omega, 
        fm = 'ml', 
        nfactors = 3)
  
  ## plots of the loadings

  corr_tools$loadings_plots <- corr_tools$fa_obj %>%  
    map(plot, 
        type = 'loadings', 
        cust_theme = globals$common_theme)
    
  corr_tools$loadings_plots <- 
    list(x = corr_tools$loadings_plots, 
         y = c('Training: IBK', 'Test: LZ/W'), 
         z = corr_tools$omega %>% 
           map_dbl(~.x$omega.tot) %>% 
           signif(2) %>% 
           paste('\u03C9 =', .)) %>% 
    pmap(function(x, y, z) x + 
           labs(title = y, 
                subtitle = z))
  
# END ----
  
  insert_tail()