# Clustering the joint (IBK + LZ) data set as suggested by the reviewer 1

  insert_head()
  
# container -------
  
  clust_rev <- list()
  
# Analysis globals -------
  
  insert_msg('Analysis globals')

  clust_rev$variables <- unique(multi_plots$coefs$variable)
  
  clust_rev$distances <- c('euclidean', 'manhattan', 'cosine')
  
  ## analysis tables: factors as numeric, 
  ## median-centered normalization
  
  clust_rev$analysis_tbl <- pah_study[c("IBK_0", "LZ_0")] %>% 
    map_dfr(select, ID, all_of(clust_rev$variables)) %>% 
    map_dfc(function(var) if(is.factor(var)) as.numeric(var) else var) %>% 
    map_dfc(function(var) if(is.numeric(var)) scale(var, center = median(var), scale = TRUE)[, 1] else var) %>% 
    column_to_rownames('ID')
  
# Clustering objects --------
  
  insert_msg('Clustering objects')
  
  clust_rev$objects[c('hcl_euclidean', 
                      'hcl_manhattan', 
                      'hcl_cosine')] <- 
    list(distance_method = clust_rev$distances, 
         k = c(2, 2, 3)) %>% 
    pmap(hcluster, 
         data = clust_rev$analysis_tbl, 
         hc_method = 'ward.D2', 
         seed = 1234)  
  
  clust_rev$objects[c('kmeans_euclidean', 
                      'kmeans_manhattan', 
                      'kmeans_cosine')] <- 
    list(distance_method = clust_rev$distances, 
         k = c(2, 2, 2)) %>% 
    pmap(kcluster, 
         data = clust_rev$analysis_tbl, 
         clust_fun = 'kmeans', 
         seed = 1234)
  
  clust_rev$objects[c('pam_euclidean', 
                      'pam_manhattan', 
                      'pam_cosine')] <- 
    list(distance_method = clust_rev$distances, 
         k = c(2, 3, 2)) %>% 
    pmap(kcluster, 
         data = clust_rev$analysis_tbl, 
         clust_fun = 'pam', 
         seed = 1234)
  
  clust_rev$objects[c('combi_euclidean', 
                      'combi_manhattan', 
                      'combi_cosine')] <- 
    list(distance_som = clust_rev$distances, 
         distance_nodes = clust_rev$distances, 
         k = c(2, 3, 3)) %>% 
    pmap(combi_cluster, 
         data = clust_rev$analysis_tbl, 
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
  
  clust_rev$variance <- clust_rev$objects %>% 
    map(var) %>% 
    map(~.x[c('total_wss', 'total_ss', 'between_ss', 'frac_var')]) %>% 
    map_dfr(as_tibble) %>% 
    mutate(algorithm = names(clust_rev$objects))
  
# Cross validation error -----
  
  insert_msg('Cross validation errors')
  
  clust_rev$cv_results <- clust_rev$objects %>% 
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
  
  clust_rev$test_results <- left_join(clust_rev$variance, 
                                      clust_rev$cv_results, 
                                      by = 'algorithm') %>% 
    left_join(., 
              clust_rev$objects %>% 
                map(extract, 'assignment') %>% 
                map_dbl(~length(unique(.x$clust_id))) %>% 
                compress(names_to = 'algorithm', 
                         values_to = 'k'), 
              by = 'algorithm') %>% 
    mutate(mean_correct = 1 - mean_error)
  
  ## plotting
  
  clust_rev$test_plot <- 
    clust_rev$test_results[c('algorithm', 'mean_correct', 'frac_var', 'k')] %>% 
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
  
# END ------
  
  insert_tail()