# Searches of the optimal clustering procedure in the exploration IBK cohort
# The cslustering features identified by elastic net regression

  insert_head()
  
# container list ----
  
  clust_dev <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals: analysis tables ------
  
  insert_msg('Globals setup')
  
  clust_dev$variables <- unique(multi_plots$coefs$variable)
  
  clust_dev$distances <- c('euclidean', 'manhattan', 'cosine')
  
  ## analysis tables: factors as numeric, 
  ## median-centered normalization
  
  clust_dev$analysis_tbl <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(select, ID, all_of(clust_dev$variables)) %>% 
    map(~map_dfc(.x, function(var) if(is.factor(var)) as.numeric(var) else var)) %>% 
    map(~map_dfc(.x, function(var) if(is.numeric(var)) scale(var, center = median(var), scale = TRUE)[, 1] else var)) %>% 
    map(column_to_rownames, 'ID')
  
# UMAP and clustering tendency ------
  
  insert_msg('UMAP and clustering tendency')
  
  ## in order to assess spontaneous clustering tendency of the dataset
  
  ## UMAP
  
  clust_dev$umap_obj <- clust_dev$analysis_tbl %>% 
    map(function(data) clust_dev$distances %>% 
          set_names(clust_dev$distances) %>% 
          map(~reduce_data(data, 
                           distance_method = .x, 
                           kdim = 2, 
                           red_fun = 'umap')))
  
  ## Hopkins stat
  
  clust_dev$clust_tendency <- 
    map2(clust_dev$analysis_tbl, 
         clust_dev$umap_obj, 
         ~c(list(data = .x), 
            map(.y, extract, 'scores') %>% 
              map(select, -observation))) %>% 
    map(~future_map(.x, 
                    get_clust_tendency, 
                    n = 50, 
                    .options = furrr_options(seed = TRUE)))

# clustering objects to be tested -----
  
  insert_msg('Clustering objects')
  
  clust_dev$objects[c('hcl_euclidean', 
                      'hcl_manhattan', 
                      'hcl_cosine')] <- 
    list(distance_method = clust_dev$distances, 
         k = c(3, 2, 2)) %>% 
    pmap(hcluster, 
         data = clust_dev$analysis_tbl$IBK_0, 
         hc_method = 'ward.D2', 
         seed = 1234)  
  
  clust_dev$objects[c('kmeans_euclidean', 
                      'kmeans_manhattan', 
                      'kmeans_cosine')] <- 
    list(distance_method = clust_dev$distances, 
         k = c(3, 3, 2)) %>% 
    pmap(kcluster, 
         data = clust_dev$analysis_tbl$IBK_0, 
         clust_fun = 'kmeans', 
         seed = 1234)
  
  clust_dev$objects[c('pam_euclidean', 
                      'pam_manhattan', 
                      'pam_cosine')] <- 
    list(distance_method = clust_dev$distances, 
         k = c(3, 3, 2)) %>% 
    pmap(kcluster, 
         data = clust_dev$analysis_tbl$IBK_0, 
         clust_fun = 'pam', 
         seed = 1234)
  
  clust_dev$objects[c('combi_euclidean', 
                      'combi_manhattan', 
                      'combi_cosine')] <- 
    list(distance_som = clust_dev$distances, 
         distance_nodes = clust_dev$distances, 
         k = c(2, 3, 3)) %>% 
    pmap(combi_cluster, 
         data = clust_dev$analysis_tbl$IBK_0, 
         xdim = 7, 
         ydim = 7,
         topo = 'hexagonal', 
         neighbourhood.fct = 'gaussian', 
         toroidal = FALSE, 
         rlen = 2000, 
         node_clust_fun = hcluster, 
         hc_method = 'ward.D2')
  
# Comparing the variances between the objects ------
  
  insert_msg('Variances for particular clustering algorithms')
  
  clust_dev$variance <- clust_dev$objects %>% 
    map(var) %>% 
    map(~.x[c('total_wss', 'total_ss', 'between_ss', 'frac_var')]) %>% 
    map_dfr(as_tibble) %>% 
    mutate(algorithm = names(clust_dev$objects))
  
# Cross validation error -----
  
  insert_msg('Cross validation errors')
  
  clust_dev$cv_results <- clust_dev$objects %>% 
    map(cv, 
        nfolds = 10, 
        kNN = 7, 
        simple_vote = FALSE, 
        resolve_ties = FALSE, 
        seed = 1234, 
        .parallel = TRUE) %>% 
    map(~.x$summary) %>% 
    compress(names_to = 'algorithm')

# A common table with the testing results, plotting -----
  
  insert_msg('Common testing result table')
  
  ## common result table, appending with the cluster number
  
  clust_dev$test_results <- left_join(clust_dev$variance, 
                                      clust_dev$cv_results, 
                                      by = 'algorithm') %>% 
    left_join(., 
              clust_dev$objects %>% 
                map(extract, 'assignment') %>% 
                map_dbl(~length(unique(.x$clust_id))) %>% 
                compress(names_to = 'algorithm', 
                         values_to = 'k'), 
              by = 'algorithm') %>% 
    mutate(mean_correct = 1 - mean_error)
  
  ## plotting
  
  clust_dev$test_plot <- 
    clust_dev$test_results[c('algorithm', 'mean_correct', 'frac_var', 'k')] %>% 
    pivot_longer(cols = all_of(c('mean_correct', 'frac_var')), 
                 names_to = 'metric',
                 values_to = 'value') %>% 
    mutate(algorithm = stri_replace(algorithm, 
                                    fixed = 'combi', 
                                    replacement = 'SOM/HCl'), 
           algorithm = stri_replace(algorithm, 
                                    fixed = 'hcl', 
                                    replacement = 'HCl'), 
           algorithm = stri_replace(algorithm, 
                                    fixed = 'kmeans', 
                                    replacement = 'k-MEANS'), 
           algorithm = stri_replace(algorithm, 
                                    fixed = 'pam', 
                                    replacement = 'PAM'), 
           algorithm = paste(algorithm, k, sep = ', k = ')) %>% 
    ggplot(aes(x = value, 
               y = reorder(stri_replace(algorithm, 
                                        fixed = '_', 
                                        replacement = ', '), 
                           value), 
               fill = metric)) + 
    geom_bar(stat = 'identity', 
             position = position_dodge(width = 0.9), 
             color = 'black') + 
    scale_fill_manual(values = c(frac_var = 'cornflowerblue', 
                                 mean_correct = 'coral3'), 
                      labels = c(frac_var = 'Frac. explained variance', 
                                 mean_correct = 'CV correct rate'), 
                      name = '') + 
    globals$common_theme + 
    theme(legend.position = 'bottom', 
          axis.title.y = element_blank()) + 
    labs(title = 'Performance of clustering algorithms', 
         subtitle = 'Clustering variance and 10-fold cross-validation')
  
# Stability of the optimal clustering solution as a function of k ------
  
  insert_msg('Stability and k')
  
  clust_dev$k_objects <- 
    list(k = 1:10) %>% 
    pmap(kcluster, 
         data = clust_dev$analysis_tbl$IBK_0, 
         distance_method = 'cosine', 
         clust_fun = 'pam', 
         seed = 1234) %>% 
    set_names(paste0('k_', 1:10))
  
  ## variance
  
  clust_dev$k_variance <- clust_dev$k_objects %>% 
    map(var) %>% 
    map_dbl(~.x$frac_var) %>% 
    compress(names_to = 'k', 
             values_to = 'frac_var') %>% 
    mutate(k = stri_extract(k, regex = '\\d{1,2}'), 
           k = as.numeric(k))
  
  ## CV stability
  
  clust_dev$k_cv <- clust_dev$k_objects %>% 
    map(cv, 
        nfolds = 10, 
        kNN = 7, 
        simple_vote = FALSE, 
        resolve_ties = FALSE, 
        seed = 1234, 
        .parallel = TRUE) %>% 
    map(~.x$summary) %>% 
    compress(names_to = 'k') %>% 
    mutate(k = stri_extract(k, regex = '\\d{1,2}'), 
           k = as.numeric(k), 
           mean_correct = 1 - mean_error)
  
  ## plotting the results
  
  clust_dev$k_plots <- 
    list(x = clust_dev[c("k_variance", "k_cv")], 
         y = c('frac_var', 'mean_correct'), 
         z = c('Clustering variance', 'CV stability'), 
         v = c('Fraction explained variance', 
               'CV correct rate')) %>% 
    pmap(function(x, y, z, v) x %>% 
           ggplot(aes(x = k, 
                      y = .data[[y]])) + 
           geom_path(color = 'steelblue') + 
           geom_point(shape = 16, 
                      size = 1.5, 
                      color = 'steelblue') + 
           geom_vline(xintercept = 2, 
                      linetype = 'dashed', 
                      color = 'coral3') +
           scale_x_continuous(breaks = 1:10) +
           globals$common_theme + 
           labs(title = z, 
                subtitle = 'PAM clustering, cosine distance', 
                x = 'Number of clusters k', 
                y = v))

# END -----
  
  plan('sequential')
  
  insert_tail()