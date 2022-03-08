# Clustering of the participants with the optimal clustering algorithm 
# chosen by cross-validation and variance fraction comparison

  insert_head()
  
# container list -----
  
  clust <- list()
  
# globals -----
  
  insert_msg('Globals setup')
  
  clust$variables <- clust_dev$variables
  
  clust$clust_obj_train <- clust_dev$objects$pam_cosine
  
  clust$clust_obj_train$clust_assignment <- extract(clust$clust_obj_train, 'assignment') %>% 
    mutate(clust_id = car::recode(clust_id, 
                                  "'1' = '#1'; '2' = '#2'; '3' = '#3'; '4' = '#4'"))
  
  clust$analysis_tbl <- clust_dev$analysis_tbl

# Prediction of the cluster assignment for the test IBK collective ----
  
  insert_msg('Prediction for the test collective')
  
  clust$clust_obj_test <- predict(clust$clust_obj_train, 
                                  newdata = clust$analysis_tbl$LZ_0, 
                                  type = 'propagation')
  
# Characteristic of the clustering objects: variance, visual diagnostic and factor importance ------
  
  insert_msg('Characteristic of the clustering objects')
  
  ## clustering variance
  
  clust$variance <- clust[c('clust_obj_train', 
                            'clust_obj_test')] %>% 
    map(var) %>% 
    map(~.x[c('total_wss', 'total_ss', 'between_ss', 'frac_var')]) %>% 
    map_dfr(as_tibble) %>% 
    mutate(cohort = c('IBK_0', 'LZ_0'))
  
  ## diagnostic plots
  
  clust$diagn_plots <- plot(clust$clust_obj_train, 
                            type = 'diagnostic', 
                            cust_theme = globals$common_theme)
  
  ## PCA plots
  
  clust$pca_plots <- clust[c('clust_obj_train', 
                             'clust_obj_test')] %>% 
    map(plot, 
        type = 'components', 
        red_fun = 'pca', 
        with = 'data', 
        k = 3, 
        cust_theme = globals$common_theme) %>% 
  map(~.x + 
        scale_fill_manual(values = globals$cluster_colors, 
                          name = 'Cluster')) %>% 
    map2(., c('Training: IBK', 'Test: LZ/W'), 
         ~.x + labs(title = .y))
  
  ## importance
  
  clust$importance <- impact(clust$clust_obj_train, 
                             seed = 1234, 
                             .parallel = TRUE)
  
# Heat maps of the clustering features -----
  
  insert_msg('Heat map of the clustering features')
  
  clust$heat_maps <- list(sample_clust_object = clust[c('clust_obj_train', 
                                                        'clust_obj_test')], 
                          plot_title = c('Training: IBK', 'Test: LZ/W')) %>% 
    pmap(plot_clust_hm, 
         cust_theme = globals$common_theme, 
         x_lab = 'Participant') %>% 
    map(~.x + 
          scale_fill_gradient2(low = 'steelblue1', 
                               mid = 'black', 
                               high = 'firebrick1', 
                               midpoint = 0, 
                               name = 'Median norm.\nfeature value', 
                               limits = c(-3, 3), 
                               oob = scales::squish) + 
          scale_y_discrete(labels = translate_vars(clust$variables), 
                           limits = rev(c('NTproBNP_log', 
                                          'PVR', 
                                          'RA_area', 
                                          'RDW_log', 
                                          'cardiac_index', 
                                          'SMWD', 
                                          'age_fc'))))
  
# Plotting of the clustering features as a panel of violin plots -----
  
  insert_msg('Violing plots of the clustering features')

  clust$vio_panels <- list(clust_object = clust[c('clust_obj_train', 
                                                'clust_obj_test')], 
                           plot_title = c('Training: IBK', 'Test: LZ/W')) %>% 
    pmap(plot_vio_panel) %>%
    map(~.x + 
          facet_grid(feature ~ ., 
                     scales = 'free', 
                     space = 'free'))

# Plotting the variable impact ------
  
  insert_msg('Plots of variable impact')
  
  clust$impact_plot <- clust$importance$summary %>% 
    ggplot(aes(x = frac_diff, 
               y = reorder(variable, frac_diff))) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed') + 
    geom_point(size = 2, 
               shape = 16, 
               color = globals$center_colors['IBK_0']) + 
    scale_y_discrete(labels = translate_vars(clust$variables)) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) +
    labs(title = 'Variable importance', 
         x = expression(Delta*' explained varaince'))
  
# END ----
  
  insert_tail()