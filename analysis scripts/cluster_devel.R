# Searches of the optimal clustering procedure in the exploration IBK cohort
# The cslustering features identified by elastic net regression

  insert_head()
  
# container list ----
  
  clust_dev <- list()
  
# globals: analysis tables ------
  
  insert_msg('Globals setup')
  
  clust_dev$variables <- unique(multi_plots$lasso_coefs$variable)
  
  clust_dev$distances <- c('euclidean', 'manhattan', 'cosine')
  
  ## analysis tables: factors as numerics, 
  ## min/max normalization
  
  clust_dev$analysis_tbl <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(select, ID, all_of(clust_dev$variables)) %>% 
    map(~map_dfc(.x, function(var) if(is.factor(var)) as.numeric(var) else var)) %>% 
    map(~map_dfc(.x, function(var) if(is.numeric(var)) scale(var, center = median(var), scale = TRUE)[, 1] else var)) %>% 
    map(column_to_rownames, 'ID')

# clustering objects to be tested -----
  
  insert_msg('Clustering objects')
  
  clust_dev$objects[c('hcl_euclidean', 
                      'hcl_manhattan', 
                      'hcl_cosine')] <- list(distance_method = clust_dev$distances, 
                                             k = c(3, 2, 2)) %>% 
    pmap(hcluster, 
         data = clust_dev$analysis_tbl$IBK_0, 
         hc_method = 'ward.D2', 
         seed = 1234)  
  
  clust_dev$objects[c('kmeans_euclidean', 
                      'kmeans_manhattan', 
                      'kmeans_cosine')] <- list(distance_method = clust_dev$distances, 
                                                k = c(3, 3, 2)) %>% 
    pmap(kcluster, 
         data = clust_dev$analysis_tbl$IBK_0, 
         clust_fun = 'kmeans', 
         seed = 1234)
  
  clust_dev$objects[c('pam_euclidean', 
                      'pam_manhattan', 
                      'pam_cosine')] <- list(distance_method = clust_dev$distances, 
                                             k = c(3, 3, 2)) %>% 
    pmap(kcluster, 
         data = clust_dev$analysis_tbl$IBK_0, 
         clust_fun = 'pam', 
         seed = 1234)
  
  clust_dev$objects[c('combi_euclidean', 
                      'combi_manhattan', 
                      'combi_cosine')] <- list(distance_som = clust_dev$distances, 
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
        nearest_n = 5, 
        seed = 1234, 
        .parallel = TRUE) %>% 
    map_dfr(~.x$summary) %>% 
    mutate(algorithm = names(clust_dev$objects))

# A common table with the testing results, plotting -----
  
  insert_msg('Common testing result table')
  
  clust_dev$test_results <- left_join(clust_dev$variance, 
                                      clust_dev$cv_results, 
                                      by = 'algorithm') %>% 
    mutate(mean_correct = 1 - mean_error)
  
  ## plotting
  
  clust_dev$test_plot <- clust_dev$test_results[c('algorithm', 'mean_correct', 'frac_var')] %>% 
    gather(key = 'metric', 
           value = 'value', 
           mean_correct, 
           frac_var) %>% 
    mutate(algorithm = stri_replace(algorithm, fixed = 'combi', replacement = 'SOM/HCl'), 
           algorithm = stri_replace(algorithm, fixed = 'hcl', replacement = 'HCl'), 
           algorithm = stri_replace(algorithm, fixed = 'kmeans', replacement = 'k-means'), 
           algorithm = stri_replace(algorithm, fixed = 'pam', replacement = 'PAM')) %>% 
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
    labs(title = 'Perfromance of clustering algorithms', 
         subtitle = 'Clustering variance and 10-fold cross-validation')
  
# END -----
  
  insert_tail()