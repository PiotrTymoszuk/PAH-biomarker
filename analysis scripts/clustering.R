# Clustering of the participants with the optimal clustering algorithm 
# chosen by cross-validation and variance fraction comparison

  insert_head()
  
# container list -----
  
  clust <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals -----
  
  insert_msg('Globals setup')
  
  clust$variables <- clust_dev$variables
  
  clust$clust_obj$IBK_0 <- clust_dev$objects$pam_cosine
  
  clust$clust_obj$IBK_0$clust_assignment <- clust$clust_obj$IBK_0 %>% 
    extract('assignment') %>% 
    mutate(clust_id = paste0('#', as.character(clust_id)), 
           clust_id = factor(clust_id))
  
  clust$analysis_tbl <- clust_dev$analysis_tbl
  
# Prediction of the cluster assignment for the test IBK collective ----
  
  insert_msg('Prediction for the test collective')
  
  clust$clust_obj$LZ_0 <- predict(clust$clust_obj$IBK_0, 
                                  newdata = clust$analysis_tbl$LZ_0, 
                                  type = 'propagation', 
                                  kNN = 7, 
                                  simple_vote = FALSE, 
                                  resolve_ties = FALSE)
  
# Characteristic of the clustering objects: variance, visual diagnostic ------
  
  insert_msg('Characteristic of the clustering objects')
  
  ## numbers of observations in the clusters
  
  clust$n_labs <- clust$clust_obj %>% 
    map(ngroups) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ', n = '))
  
  ## clustering variance
  
  clust$variance <- clust$clust_obj %>% 
    map(var) %>% 
    map(~.x[c('total_wss', 'total_ss', 'between_ss', 'frac_var')]) %>% 
    map(as_tibble) %>% 
    compress(names_to = 'cohort')
  
  ## diagnostic plots
  
  clust$diagn_plots <- plot(clust$clust_obj$IBK_0, 
                            type = 'diagnostic', 
                            cust_theme = globals$common_theme)
  
  ## distance heat maps
  
  clust$diagn_plots$dist_hm <- clust$clust_obj %>% 
    map(plot, 
        type = 'heat_map', 
        cust_theme = globals$common_theme) %>% 
    map2(., c('Training: IBK', 'Test: LZ/W'), 
         ~.x + 
           labs(title = .y) + 
           theme(axis.text = element_blank(),
                 axis.text.x = element_blank(), 
                 axis.line = element_blank(), 
                 axis.title = element_blank()))
  
  ## distance MDS
  
  clust$diagn_plots$dist_mds <- clust$clust_obj %>% 
    map(plot, 
        type = 'components',
        red_fun = 'mds', 
        kdim = 2, 
        cust_theme = globals$common_theme) %>% 
    map2(., clust$n_labs, 
         ~.x + 
          scale_fill_manual(values = globals$cluster_colors, 
                            labels = .y, 
                            name = 'Cluster')) %>% 
    map2(., c('Training: IBK', 'Test: LZ/W'), 
         ~.x + labs(title = .y))
  
# PCA and UMAP projections --------
  
  insert_msg('PCA and UMAP plots')
  
  ## PCA plots
  
  clust$pca_plots <- clust$clust_obj %>% 
    map(plot, 
        type = 'components', 
        red_fun = 'pca', 
        with = 'data', 
        k = 2, 
        cust_theme = globals$common_theme) %>% 
    map2(., clust$n_labs, 
         ~.x + 
           scale_fill_manual(values = globals$cluster_colors, 
                             labels = .y, 
                             name = 'Cluster')) %>% 
    map2(., c('Training: IBK', 'Test: LZ/W'), 
         ~.x + labs(title = .y))
  
  ## UMAP plots
  
  clust$umap_plots <- clust$clust_obj %>% 
    map(plot, 
        type = 'components', 
        red_fun = 'umap', 
        with = 'data', 
        k = 2, 
        cust_theme = globals$common_theme, 
        random_state = 123) %>% 
    map2(., clust$n_labs, 
         ~.x + 
           scale_fill_manual(values = globals$cluster_colors, 
                             labels = .y, 
                             name = 'Cluster')) %>% 
    map2(., c('Training: IBK', 'Test: LZ/W'), 
         ~.x + labs(title = .y))
  
# Heat maps of the clustering features -----
  
  insert_msg('Heat map of the clustering features')
  
  clust$heat_maps <- 
    list(x_object = clust$clust_obj, 
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
          scale_y_discrete(labels = exchange(clust$variables, 
                                             dict = globals$var_labs), 
                           limits = rev(c('NTproBNP_log', 
                                          'PVR', 
                                          'RA_area', 
                                          'RDW_log', 
                                          'cardiac_index', 
                                          'SMWD', 
                                          'age_fc'))))
  
# Plotting of the clustering features as a panel of violin plots -----
  
  insert_msg('Violin plots of the clustering features')

 clust$vio_panels <- 
   list(clust_object = clust$clust_obj, 
        plot_title = c('Training: IBK', 'Test: LZ/W')) %>% 
   pmap(plot_vio_panel) %>%
   map(~.x + 
         facet_grid(feature ~ ., 
                    scales = 'free', 
                    space = 'free'))

# Variable importance ------
  
  insert_msg('Variable importance')
  
  ## permutation importance
 
  set.seed(1234)

  clust$importance <- sample(1:1000, 100, replace = FALSE) %>% 
    set_names(paste0('run_', 1:100)) %>% 
    future_map(~impact(clust$clust_obj$IBK_0, 
                       seed = .x, 
                       .parallel = FALSE), 
               .options = furrr_options(seed = TRUE))

  
  ## formatting the results
  
  clust$importance <- clust$importance %>% 
    map(select, variable, frac_diff) %>% 
    compress(names_to = 'run') %>% 
    filter(variable != 'data')
  
# Plotting the variable impact ------
  
  insert_msg('Plots of variable impact')
  
  clust$impact_plot <- clust$importance %>% 
    ggplot(aes(x = frac_diff, 
               y = reorder(variable, frac_diff))) + 
    geom_violin(fill = globals$center_colors["IBK_0"], 
                alpha = 0.25) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed') + 
    geom_point(size = 2, 
               shape = 16, 
               alpha = 0.45, 
               color = globals$center_colors['IBK_0'], 
               position = position_jitter(width = 0, 
                                          height = 0.1)) + 
    scale_y_discrete(labels = exchange(clust$variables, 
                                       dict = globals$var_labs)) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) +
    labs(title = 'Variable importance', 
         subtitle = 'n = 100 random re-shuffles, IBK cohort', 
         x = expression(Delta*' explained variance'))
  
# END ----
  
  plan('sequential')
  
  insert_tail()