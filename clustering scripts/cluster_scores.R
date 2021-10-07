# This script compares the risk scores (developed and comparators) between the clinical clusters

  insert_head()
  
# data container -----
  
  clust_scores <- list()
  
# globals -----
  
  insert_msg('Globals setup')
  
  clust_scores$variables <- c(multi_modeling$cv_pass$stats$model_id, 
                              pah_study$comparators$variable)
  
# serial analysis -----
  
  insert_msg('Serial analysis')
  
  ## analysis objects
  
  clust_scores$analyses <- study_clust$analysis_tbl %>% 
    map(function(cohort) clust_scores$variables %>% 
          map(safely(analyze_feature), 
              inp_tbl = cohort, 
              split_var = 'clust_id') %>% 
          set_names(clust_scores$variables) %>% 
          transpose %>% 
          map(compact)) %>% 
    map(~.x$result) %>% 
    set_names(names(study_clust$analysis_tbl))
  
  ## summaries
  
  clust_scores$summaries <- clust_scores$analyses %>% 
    map(function(cohort) cohort %>% 
          map_dfr(extract_test_summary) %>% 
          filter(test == 'kruskal') %>% 
          mutate(p_adj = p.adjust(p_value, 'BH'), 
                 significance = ifelse(p_adj >= 0.05,
                                       'ns',
                                       paste(signif(p_adj, 2)))))
  
# Plotting: violinos, appending with the adjusted p values -----
  
  insert_msg('Plotting the results')
  
  ## base plots
  
  clust_scores$plots <-  c('IBK_0', 'LZ_0') %>% 
    map(function(x) list(analysis_obj = clust_scores$analyses[[x]], 
                         label = ifelse(names(clust_scores$analyses[[x]]) %in% pah_study$comparators$variable, 
                                        globals$comp_labs[names(clust_scores$analyses[[x]])], 
                                        stri_replace(names(clust_scores$analyses[[x]]), fixed = 'sign_', 'Signature ')) %>% 
                           paste(., globals$center_labs[[x]], sep = ': ')) %>% 
          pmap(plot_analysis, 
               violin = T, 
               labeller = study_clust$clust_labs, 
               fill_colors = unname(study_clust$clust_colors), 
               cust_theme = globals$common_theme, 
               y_lab = 'Risk score', 
               point_alpha = 0.8) %>% 
          set_names(names(clust_scores$analyses[[x]]))) %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
  ## p values
  
  clust_scores$plots <- c('IBK_0', 'LZ_0') %>% 
    map(function(x) list(plot = clust_scores$plots[[x]], 
                         sign = clust_scores$summaries[[x]][['significance']]) %>% 
          pmap(function(plot, sign) plot + 
                 labs(subtitle = sign))) %>% 
    set_names(c('IBK_0', 'LZ_0'))
  
# Plotting table with the median min/max normalized score values for the signature ------
  
  insert_msg('Table with median score values in the clusters')

  clust_scores$median_tbl <- study_clust$analysis_tbl %>% 
    map(select, 
        all_of(pah_study$comparators$variable), 
        all_of(c('sign_309', 'sign_2525'))) %>% 
    map(function(cohort) cohort %>% 
          map_dfc(min_max)) %>% 
    map2(map(study_clust$analysis_tbl , ~.x[, 'clust_id']), ., cbind)
  
  clust_scores$median_tbl <- clust_scores$median_tbl %>% 
    map(gather, 
        key = 'variable', 
        value = 'score', 
        all_of(pah_study$comparators$variable), 
        all_of(c('sign_309', 'sign_2525'))) %>% 
    map(group_by, 
        variable, 
        clust_id) %>% 
    map(summarise, 
        score = median(score, na.rm = T))
  
  clust_scores$median_tbl <- map2(clust_scores$median_tbl, 
                                  map(clust_scores$summaries, 
                                      ~.x[c('variable', 'p_adj')]), 
                                  left_join, 
                                  by = 'variable')
  
  clust_scores$median_tbl <- clust_scores$median_tbl %>% 
    map2_dfr(., names(.), ~mutate(.x, cohort = .y)) %>% 
    mutate(clust_name = study_clust$clust_labs[clust_id], 
           var_label = ifelse(variable %in% pah_study$comparators$variable, 
                              globals$comp_labs[variable], 
                              stri_replace(variable, fixed = 'sign_', replacement = 'Sign ')), 
           var_type = ifelse(variable %in% pah_study$comparators$variable, 
                             'comparator', 'new') %>% 
             factor(c('new', 'comparator')))
  
# Representation of the median score in the cluster as a heat map -----
  
  insert_msg('Heat map of the ,edian scores')

  clust_scores$heat_map <- clust_scores$median_tbl %>% 
    ggplot(aes(x = clust_name, 
               y = reorder(var_label, -p_adj), 
               fill = score)) + 
    geom_tile(color = 'gray60') + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'white', 
                         high = 'firebrick', 
                         midpoint = 0.5, 
                         limits = c(0, 1), 
                         name = 'Normalized\nmedian') + 
    facet_grid(var_type ~ cohort, 
               labeller = labeller(.cols = globals$center_labs, 
                                   .rows = globals$signature_labels), 
               scales = 'free_y', 
               space = 'free') + 
    globals$common_theme + 
    theme(axis.title = element_blank(), 
          axis.line = element_blank()) + 
    labs(title = 'Risk score values in participant clusters', 
         subtitle = 'Min/max-normaized median score', 
         tag = study_clust$plot_tags %>% 
           map2_chr(globals$center_labs[c('IBK_0', 'LZ_0')], ., 
                    paste0, 
                    collapse = ': '))
  
# Representation of the significance as a bar plot -----
  
  insert_msg('Bar plot with significances')
  
  clust_scores$signif_bar <- clust_scores$median_tbl %>% 
    filter(clust_id == 'cluster_1') %>% 
    mutate(sign = ifelse(p_adj < 0.05, 'yes', 'no')) %>% 
    ggplot(aes(x = -log10(p_adj), 
               y = reorder(var_label, -p_adj), 
               fill = sign))  + 
    geom_bar(stat = 'identity', 
             color = 'gray60') + 
    geom_vline(xintercept = -log10(0.05), 
               linetype = 'dashed') + 
    scale_fill_manual(values = c(no = 'gray60', 
                                 yes = 'coral3')) + 
    guides(fill = F) +
    facet_grid(var_type ~ cohort, 
               labeller = labeller(.cols = globals$center_labs, 
                                   .rows = globals$signature_labels), 
               scales = 'free_y', 
               space = 'free') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank(), 
          panel.grid.major = element_line(color = 'gray90')) + 
    labs(title = 'Risk score values in participant clusters', 
         subtitle = 'Kruskal-Wallis test', 
         x = expression('-log'[10]*'pFDR'))
  
# END -----
  
  insert_tail()