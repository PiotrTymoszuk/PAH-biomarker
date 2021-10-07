## a script providing tools for hierarchical and kmeans or PAM clustering

# tools -----

  require(tidyverse)
  require(philentropy)
  require(factoextra)
  require(dendextend)
  require(scrime)
  require(cluster)
  require(pcaPP)
  require(dbscan)

# helper functions -----

  calculate_dist <- function(inp_tbl, method) {
    
    ## calculates a dissimilarity matrix
    
    if(method == 'smc') {
      
      dist_mtx <- inp_tbl %>% 
        as.matrix %>% 
        smc(dist = T)
      
      
    } else {
      
      dist_mtx <- inp_tbl %>% 
        as.matrix %>% 
        distance(method = method, 
                 use.row.names = T)
      
    }
    
    return(dist_mtx)
    
  }
  
  plot_knn_distance <- function(diss_obj, k, eps = NULL) {
    
    ## makes k-NN distance plot
    
    sort_distances <- kNNdist(x = diss_obj, 
                              k = k) %>% 
      sort %>% 
      tibble(sample = names(.), 
             knn_dist = unname(.)) %>% 
      select(sample, 
             knn_dist) %>% 
      mutate(point = 1:nrow(.))
    
    dist_plot <- sort_distances %>% 
      ggplot(aes(x = point, 
                 y = knn_dist)) + 
      geom_area(color = 'cornsilk4', 
                fill = 'cornsilk2') + 
      theme_classic() + 
      labs(x = 'Sample', 
           y = paste(k, 'NN distance', sep = '-'))
    
    if(!is.null(eps)) {
      
      dist_plot <- dist_plot + 
        geom_hline(yintercept = eps, 
                   linetype = 'dashed', 
                   color = 'black')
      
    }
    
    return(dist_plot)
    
  }

# dendrogram plotting -----

  plot_dendro <- function(clust_str, k, labels = T, cluster_colors = NULL, cluster_labels = NULL, 
                          cluster_leg_title = NULL, plot_title = NULL, plot_subtitle = NULL, y_lab = NULL, ...) {
    
    ## a wrapper for the respective function from the factoextra package
    
    dendro_plot <- as.dendrogram(clust_str) %>% 
      color_branches(k = k) %>% 
      set('branches_lwd', 0.5) %>% 
      set('labels_cex', 0.5) %>%
      ggplot(labels = labels) + 
      theme_classic() + 
      theme(axis.line.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.title.x = element_blank()) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           y = y_lab)
    
    if(all(c(!is.null(cluster_colors), 
             !is.null(cluster_labels), 
             !is.null(cluster_leg_title)))) {
      
      suppressMessages(dendro_plot <- dendro_plot + 
                         scale_color_manual(values = cluster_colors, 
                                            labels = cluster_labels, 
                                            name = cluster_leg_title) + 
                         guides(color = guide_legend()))
      
    }
    
    return(dendro_plot)
    
  }


# hierarchical, kmeans, pam and dbscan clustering -----

  hcluster_data <- function(inp_tbl, distance_method = 'jaccard', k = 2, hc_method = 'ward.D2', seed = 1234, ...) {
    
    ## calculates distance between the data in the data frame with named rows using philentropy package
    ## the simple matching distance id provided by scrime
    ## and subjects the distance matrix to clustering using the hierarchical clustering algorithm
    ## provides plot of mean squared error vs cluster number and silhouette stat vs cluster number 
    ## to facilitate the decision on the cluster number
    
    set.seed(seed = seed)
    
    ## distance calculation
    
    dist_mtx <- calculate_dist(inp_tbl = inp_tbl, 
                               method = distance_method)
    
    ## diagnostic plots for cluster numbers
    
    diagn_hclust_wss <- fviz_nbclust(dist_mtx, 
                                     FUNcluster = hcut,
                                     hc_method = hc_method, 
                                     method = 'wss') + 
      geom_vline(xintercept = k, 
                 linetype = 'dashed', 
                 color = 'steelblue') + 
      theme_classic() + 
      labs(title = 'Optimal cluster number', 
           subtitle = paste('Hierarchical clustering,', 
                            distance_method, 
                            'distance,', 
                            hc_method, 
                            'method'))
    
    diagn_hclust_silhouette <- fviz_nbclust(dist_mtx, 
                                            FUNcluster = hcut,
                                            hc_method = hc_method, 
                                            method = 'silhouette') + 
      theme_classic() + 
      labs(title = 'Optimal cluster number', 
           subtitle = paste('Hierarchical clustering,', 
                            distance_method, 
                            'distance,', 
                            hc_method, 
                            'method'))
    
    ## hierarchical clustering and cluster assignment table
    
    hclust_str <- hclust(dist_mtx %>% 
                           as.dist, 
                         method = hc_method, ...)
    
    hclust_ass <- hclust_str %>% 
      cutree(k = k) %>% 
      tibble(variable = names(.), 
             clust_id = unname(.)) %>% 
      select(variable, 
              clust_id)
    
    ## dendrogram plot for hierarchical clustering
    
    dendro_plot <- plot_dendro(clust_str = hclust_str, 
                               k = k, 
                               labels = T)
    
    clust_analysis <- list(dist_tbl = dist_mtx %>%  as.data.frame, 
                           dist_obj = dist_mtx %>% as.dist, 
                           dist_method = distance_method, 
                           clust_fun = 'hclust', 
                           diagnostic_plots = list(wss = diagn_hclust_wss, 
                                                   silhouette = diagn_hclust_silhouette, 
                                                   dendrogram = dendro_plot), 
                           clust_obj = hclust_str, 
                           clust_assignment = hclust_ass)
    
    attr(clust_analysis, 'class') <- 'cluster_analysis'
    
    return(clust_analysis)
    
  }
  
  kcluster_data <- function(inp_tbl, distance_method = 'jaccard', clust_fun = kmeans, k = 2, seed = 1234, ...) {
    
    ## calculates distance between the data in the data frame with named rows using philentropy package
    ## the simple matching distance is provided by scrime
    ## and subjects the distance matrix to clustering using the kmeans or pam clustering algorithm
    ## provides plot of mean squared error vs cluster number and silhouette stat vs cluster number 
    ## to facilitate the decision on the cluster number
    
    set.seed(seed = seed)
    
    ## distance calculation
    
    dist_mtx <- calculate_dist(inp_tbl = inp_tbl, 
                               method = distance_method)
    
    ## diagnostic plots for cluster numbers
    
    diagn_kclust_wss <- fviz_nbclust(dist_mtx, 
                                     FUNcluster = clust_fun, 
                                     method = 'wss') + 
      geom_vline(xintercept = k, 
                 linetype = 'dashed', 
                 color = 'steelblue') + 
      theme_classic() + 
      labs(title = 'Optimal cluster number', 
           subtitle = paste(quote(clust_fun) %>% as.character, 
                            'clustering,', 
                            distance_method, 
                            'distance'))
    
    diagn_kclust_silhouette <- fviz_nbclust(dist_mtx, 
                                            FUNcluster = clust_fun, 
                                            method = 'silhouette') + 
      theme_classic() + 
      labs(title = 'Optimal cluster number', 
           subtitle = paste(quote(clust_fun) %>% as.character, 
                            'clustering,', 
                            distance_method, 
                            'distance'))
    
    ## kmeans/pam clustering and cluster assignment table

    kclust_str <- clust_fun(dist_mtx %>% 
                           as.dist, 
                           k, ...)
    
    kclust_ass <- kclust_str$cluster %>% 
      tibble(variable = names(.), 
             clust_id = unname(.)) %>% 
      select(variable, 
             clust_id)
      
    
    ## output
    
    clust_analysis <- list(dist_tbl = dist_mtx %>% as.data.frame, 
                           dist_obj = dist_mtx %>% as.dist, 
                           dist_method = distance_method, 
                           clust_fun = quote(clust_fun) %>% as.character, 
                           diagnostic_plots = list(wss = diagn_kclust_wss, 
                                                   silhouette = diagn_kclust_silhouette), 
                           clust_obj = kclust_str, 
                           clust_assignment = kclust_ass)
    
    attr(clust_analysis, 'class') <- 'cluster_analysis'
    
    return(clust_analysis)
    
  }
  
  dbscan_data <- function(inp_tbl, distance_method = 'euclidean', eps, minPts = 5, ...) {
    
    ## calculates distance between the data in the data frame with named rows using philentropy package, 
    ## the simple matching distance is provided by scrime
    ## and subjects the distance matrix to clustering using the dbscan clustering algorithm
    ## to facilitate the decision on the optimal eps value, a plot of kNN distance (where k = minPts - 1)
    ## versus sample is provided

    ## distance calculation
    
    dist_mtx <- calculate_dist(inp_tbl = inp_tbl, 
                               method = distance_method)
    
    ## diagnostic plot
    
    knn_dist_plot <- plot_knn_distance(dist_mtx %>% 
                                         as.dist, 
                                       k = minPts - 1, 
                                       eps = eps)
    
    ## dbscan clustering and cluster assignment table
    
    kclust_str <- dbscan(dist_mtx %>% 
                              as.dist, 
                         eps = eps, 
                         minPts = minPts, ...)
    
    kclust_ass <- tibble(variable = rownames(dist_mtx), 
                         clust_id = kclust_str$cluster)
    
    
    ## output
    
    clust_analysis <- list(dist_tbl = dist_mtx %>% as.data.frame, 
                           dist_obj = dist_mtx %>% as.dist, 
                           dist_method = distance_method, 
                           clust_fun = 'dbscan', 
                           diagnostic_plots = list(knn_distance = knn_dist_plot), 
                           clust_obj = kclust_str, 
                           clust_assignment = kclust_ass)
    
    attr(clust_analysis, 'class') <- 'cluster_analysis'

    return(clust_analysis)    
    
  }

# plotting clustering assignment as a heat map of the distances -----
  
  plot_dist_hm <- function(cluster_analysis = NULL, dist_tbl = NULL, clust_assign_tbl = NULL, 
                           dist_transf_fun = function(x) x, 
                           cluster_labs = NULL, cust_theme = NULL, min_color = 'steelblue4', 
                           max_color = 'firebrick4', mid_color = 'white', mid_point = NULL, 
                           legend_title = 'distance', plot_title = NULL, plot_subtitle = NULL, plot_tag = NULL) {
    
    ## makes a heat plot plot of the distances between the observartions/variables
    ## with a faceting by cluster structure
    ## extracts the data from the clustr_analysis object or from a dissimilarity table
    
    
    
    ## plotting table
    
    if(is.null(dist_tbl) & is.null(clust_assign_tbl)) {
      
      dist_tbl <- cluster_analysis$dist_tbl
      
      clust_assign_tbl <- cluster_analysis$clust_assignment
      
    }
    
    clust_vars <- names(dist_tbl)
    
    plotting_tbl <- dist_tbl %>% 
      rownames_to_column('variable') %>%
      left_join(clust_assign_tbl, 
                by = 'variable') %>% 
      gather(key = 'variable2', 
             value = distance, 
             all_of(clust_vars)) %>% 
      mutate(distance = dist_transf_fun(distance))
    
    plotting_tbl <- left_join(plotting_tbl, 
                              set_names(clust_assign_tbl[, c(1:2)], 
                                        c('variable2', 
                                          'clust_id2')), 
                              by = 'variable2')
    
    ## plotting
    
    if(is.null(mid_point)) {
      
      midpoint <- range(plotting_tbl$distance) %>% 
        mean
      
    }
    
    if(is.null(cust_theme)) {
      
      cust_theme <- theme_minimal()
      
    }
    
    heat_map <- plotting_tbl %>% 
      ggplot(aes(x = reorder(variable, clust_id), 
                 y = reorder(variable2, clust_id2), 
                 fill = distance)) + 
      geom_tile() + 
      scale_fill_gradient2(low = min_color, 
                           mid = mid_color, 
                           high = max_color, 
                           midpoint = midpoint) + 
      cust_theme + 
      theme(axis.title = element_blank(), 
            axis.text.x = element_text(angle = 90, 
                                       hjust = 1, 
                                       vjust = 0.5)) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           fill = legend_title)
    
    if(is.null(cluster_labs)) {
      
      heat_map <- heat_map + 
        facet_grid(clust_id2 ~ clust_id, 
                   scales = 'free', 
                   space = 'free')
      
    } else {
      
      heat_map <- heat_map + 
        facet_grid(clust_id2 ~ clust_id, 
                   scales = 'free', 
                   space = 'free', 
                   labeller = labeller(.cols = cluster_labs, 
                                       .rows = cluster_labs))
      
    }
    
  }
  
  plot_clust_mds <- function(cluster_analysis, k_dim = 2, red_fun = 'mds', 
                             cluster_labs = NULL, cluster_colors = NULL, cust_theme = NULL, 
                             plot_title = NULL, plot_subtitle = NULL, plot_tag = NULL, 
                             jitter_width = 0, jitter_height = 0, legend_title = 'Cluster', ...) {
    
    ## makes k-dimension multi-dimensional scaling or PCA of the dissimilarity matrix 
    ## and plots the first two dimensions with the color coding for the clusters
    
    ## dimensionality reduction
    
    if(red_fun == 'mds') {
      
      plotting_tbl <- cluster_analysis$dist_tbl %>% 
        cmdscale(k = k_dim, ...)
      
      red_obj <- NULL
      
    } else {
      
      id_vec <- rownames(cluster_analysis$dist_tbl)
      
      red_obj <- cluster_analysis$dist_tbl %>% 
        PCAproj(k = k_dim)
      
      plotting_tbl <- red_obj$scores
      
      rownames(plotting_tbl) <- id_vec
      
    }
    
    plotting_tbl <- plotting_tbl %>% 
      as.data.frame %>% 
      set_names(paste('dim', 1:ncol(.), sep = '_')) %>% 
      rownames_to_column('variable') %>% 
      as_tibble
    
    ## appending the dimension table with the clustering information
    
    plotting_tbl <- plotting_tbl %>% 
      left_join(cluster_analysis$clust_assignment, 
                by = 'variable')
    
    ## plotting
    
    if(is.null(cust_theme)) {
      
      cust_theme <- theme_classic()
      
    }
    
    point_plot <- plotting_tbl %>% 
      ggplot(aes(x = dim_1, 
                 y = dim_2, 
                 fill = factor(clust_id))) + 
      geom_point(shape = 21, 
                 size = 2, 
                 position = position_jitter(width = jitter_width, 
                                            height = jitter_height)) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = 'Dimension 1', 
           y = 'Dimension 2', 
           fill = legend_title) + 
      cust_theme
    
    if(!is.null(cluster_labs) & !is.null(cluster_colors)) {
      
      point_plot <- point_plot  + 
        scale_fill_manual(values = cluster_colors, 
                          labels = cluster_labs, 
                          name = legend_title)
      
    }
    
    return(list(plot = point_plot, 
                red_obj = red_obj))
    
  }
  
# END ----