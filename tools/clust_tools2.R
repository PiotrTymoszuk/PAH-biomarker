## a script providing tools for hierarchical, kmeans, PAM, DB scan and SOM clustering

# tools -----

  require(tidyverse)
  require(philentropy)
  require(factoextra)
  require(dendextend)
  require(nomclust)
  require(cluster)
  require(pcaPP)
  require(dbscan)
  require(rlang)
  require(kohonen)
  require(Rcpp)
  require(grDevices)
  require(cowplot)
  require(ggrepel)
  require(coxed)
  require(caret)
  require(furrr)

  ## additional kernels for kohonen::som()

  tryCatch(library(somKernels), 
           error = function(e) install.packages('./tools/somKernels.tar', 
                                                repos = NULL, 
                                                type = 'source', 
                                                INSTALL_opts = '--no-help'))

  map <- purrr::map

  map_som <- kohonen::map
  
# helper functions -----

  calculate_dist <- function(data, method) {
    
    ## calculates a dissimilarity matrix
    
    if(!is_tibble(data) & !is.data.frame(data) & !is.matrix(data)) {
      
      stop('Provide a valid data.frame, tibble or matrix as input data', call. = FALSE)
      
    }
    
    av_distances <- c(getDistMethods(), 
                      'smc', 
                      'sumofsquares')
    
    if(method == 'smc') {
      
      dist_mtx <- data %>% 
        as.data.frame %>% 
        sm %>% 
        as.matrix
      
    } else {
      
      if(method == 'sumofsquares') method <- 'squared_euclidean'
      
      dist_mtx <- data %>% 
        as.matrix %>% 
        distance(method = method, 
                 use.row.names = TRUE, 
                 mute.message = TRUE)
      
      if(method %in% c('cosine', 'ruzicka')) {
        
        ### handling the 0 - 1 similarity coefficients
        
        dist_mtx <- 1 - dist_mtx
        
      }
      
    }
    
    return(dist_mtx)
    
  }
  
  plot_knn_distance <- function(diss_obj, 
                                k, 
                                eps = NULL, 
                                plot_title = NULL, 
                                plot_subtitle = NULL, 
                                plot_tag = NULL, 
                                cust_theme = theme_classic()) {
    
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
      cust_theme + 
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
  
  plot_dendro <- function(clust_str, 
                          k, 
                          labels = T, 
                          cluster_colors = NULL, 
                          cluster_labels = NULL, 
                          cluster_leg_title = NULL, 
                          plot_title = NULL, 
                          plot_subtitle = NULL, 
                          plot_tag = NULL, 
                          y_lab = NULL, 
                          cust_theme = theme_classic(), ...) {
    
    ## a wrapper for the respective function from the factoextra package
    
    dendro_plot <- as.dendrogram(clust_str) %>% 
      color_branches(k = k) %>% 
      set('branches_lwd', 0.5) %>% 
      set('labels_cex', 0.5) %>%
      ggplot(labels = labels) + 
      cust_theme + 
      theme(axis.line.x = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.title.x = element_blank()) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
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
  
  plot_nbclust <- function(dist_mtx, 
                           k, 
                           FUNcluster = NULL, 
                           method = 'wss', 
                           plot_title = NULL, 
                           plot_subtitle = NULL, 
                           plot_tag = NULL, 
                           cust_theme = theme_classic(), ...) {
    
    ## plots a WSS curve as a function of k
    
    nb_plot <- fviz_nbclust(dist_mtx, 
                            FUNcluster = FUNcluster, 
                            method = method, ...) + 
      geom_vline(xintercept = k, 
                 linetype = 'dashed', 
                 color = 'coral3') + 
      cust_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag)
    
    return(nb_plot)
    
  }
  
  plot_som <- function(kohonen_object) {
    
    ## generates diagnostic plots for the test cohort
    
    if(class(kohonen_object) != 'kohonen') {
      
      stop('The input has to be a kohonen object', call. = FALSE)
      
    }
    
    plot_exprs <- enexpr(kohonen_object)
    
    plot_exprs_lst <- c('change', 
                        'codes', 
                        'counts', 
                        'mapping', 
                        'dist.neighbours', 
                        'quality') %>% 
      purrr::map(function(x) expr(~plot(!!plot_exprs, !!x))) %>% 
      set_names(c('change', 
                  'codes', 
                  'counts', 
                  'mapping', 
                  'dist.neighbours', 
                  'quality'))
    
    diagn_plots <- plot_exprs_lst %>% 
      purrr::map(function(x) ggdraw(as_grob(eval_tidy(x))))
    
    return(diagn_plots)
    
  }
  
  plot_train_som <- function(kohonen_object, 
                             plot_title = NULL, 
                             plot_subtitle = NULL, 
                             cust_theme = theme_classic(), ...) {
    
    ## plots the distance to the winning unit (neuron) during the training process
    ## a nicer plot using ggplot(). May pass additional arguments to the smoothing function
    
    if(class(kohonen_object) != 'kohonen') {
      
      stop('The training plots available only for the SOM cluster analyses', 
           call. = F)
      
    }
    
    change_plot <- kohonen_object$changes %>% 
      as.data.frame %>% 
      set_names('dist') %>% 
      rownames_to_column('Iteration') %>% 
      ggplot(aes(x = as.numeric(Iteration), 
                 y = dist)) + 
      geom_point(shape = 16, 
                 size = 1, 
                 color = 'gray60') + 
      geom_smooth(color = 'steelblue', 
                  method = 'loess', ...) + 
      cust_theme + 
      labs(x = 'Iteration', 
           y = 'Mean distance\nto the winning unit', 
           title = plot_title, 
           subtitle = plot_subtitle, 
           tag = paste0('\nIterations: n = ', 
                        nrow(kohonen_object$changes)))
    
    return(change_plot)
    
  }
  
  plot_point <- function(data, 
                         x_var, 
                         y_var, 
                         fill_var = NULL, 
                         label_var = NULL, 
                         plot_title = NULL, 
                         plot_subtitle = NULL, 
                         plot_tag = NULL, 
                         x_lab = x_var, 
                         y_lab = y_var, 
                         fill_lab = NULL, 
                         cust_theme = theme_classic(), 
                         point_color = 'steelblue', 
                         point_alpha = 1, 
                         show_segments = FALSE, 
                         segment_color = 'steelblue', 
                         segment_alpha = 1, 
                         label_color = point_color, 
                         txt_color = 'black', 
                         txt_size = 2.5, 
                         txt_type = c('label', 'text'), 
                         jitter_width = 0, 
                         jitter_height = 0) {
    
    ## makes a highly customized point plot
    
    if(is.null(fill_var)) {
      
      pplot <- data %>% 
        ggplot(aes(x = .data[[x_var]], 
                   y = .data[[y_var]]))
      
    } else {

      pplot <- data %>% 
        ggplot(aes(x = .data[[x_var]], 
                   y = .data[[y_var]], 
                   fill = .data[[fill_var]]))
      
    }
    
    if(show_segments) {
      
      pplot <- pplot + 
        geom_segment(aes(x = 0, 
                         y = 0, 
                         xend = .data[[x_var]], 
                         yend = .data[[y_var]]), 
                     color = segment_color, 
                     alpha = segment_alpha)
      
    }
    
    if(is.null(fill_var)) {
      
      pplot <- pplot + 
        geom_point(shape = 21, 
                   size = 2, 
                   fill = point_color, 
                   alpha = point_alpha, 
                   position = position_jitter(width = jitter_width, 
                                              height = jitter_height))
      
    } else {
      
      pplot <- pplot + 
        geom_point(shape = 21, 
                   size = 2, 
                   alpha = point_alpha, 
                   position = position_jitter(width = jitter_width, 
                                              height = jitter_height))
      
    }
    
    if(!is.null(label_var)) {
      
      txt_type <- match.arg(txt_type[1], c('label', 'text'))
      
      if(txt_type == 'label') {
        
        pplot <- pplot + 
          geom_label_repel(aes(label = .data[[label_var]]), 
                           size = txt_size, 
                           fill = label_color, 
                           color = txt_color, 
                           box.padding = 0.1, 
                           label.padding = 0.1)
        
      } else {
        
        pplot <- pplot + 
          geom_text_repel(aes(label = .data[[label_var]]), 
                          size = txt_size, 
                          color = txt_color, 
                          box.padding = 0.1)
        
      }
      
    }
    
    pplot + 
      cust_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = x_lab, 
           y = y_lab, 
           fill = fill_lab)
      
    
  }
  
  get_sum_sq <- function(dist_mtx, assignment) {
    
    ## calculates within and total sum of squares for any clustering object
    
    ## within sum of squares
    
    x_ss <- aggregate(dist_mtx, 
                      by = assignment[, 'clust_id'], 
                      FUN = function(x) sum(scale(x, scale = FALSE)^2))
    
    wss <- rowSums(x_ss[, -1])
    
    total_wss <- sum(x_ss[, -1])
    
    ## total sum of squares
    
    total_ss <- sum(scale(dist_mtx, scale = FALSE)^2)
    
    ## output
    
    list(wss = wss, 
         total_wss = total_wss, 
         total_ss = total_ss, 
         between_ss = total_ss - total_wss, 
         frac_var = (total_ss - total_wss)/total_ss)
    
  }
  
  get_data_dim <- function(data) {
    
    list(observations = nrow, 
         variables = ncol) %>% 
      map(~.x(data))
    
  }
  
  vote_simple <- function(vector) {
    
    ## returns the majority class present in a vector
    
    voting_res <- sort(table(vector), decreasing = TRUE)
    
    names(voting_res)[1]
    
  }
  
  vote_kernel <- function(vector, dist_vec, kernel_fun = function(x) 1/x) {
    
    ## distance weighted voting
    
    vote_tbl <- tibble(raw_votes = vector, 
                       distances = dist_vec, 
                       weighted_votes = kernel_fun(dist_vec)) %>% 
      filter(complete.cases(.)) %>% 
      group_by(raw_votes) %>% 
      summarise(vote_sum = sum(weighted_votes)) %>% 
      arrange(-vote_sum)
    
    vote_tbl$raw_votes[1]
    
  }
  
# OOP definitions -----
  
  clust_analysis <- function(x) {
    
    ## creates a clust_analysis object
    ## given a result list
    
    if(!is.list(x)) stop('Only named lists can be converted to a clust_analysis objects', call. = FALSE)
    
    if(is.null(names(x))) stop('Only named lists can be converted to a clust_analysis objects', call. = FALSE)
    
    if(!all(c('data', 'dist_mtx', 'dist_method', 'clust_obj', 'clust_fun', 'clust_assignment') %in% names(x))) {
      
      stop('The input list needs to provide data, dist_mtx, dist_method, clust_fun, clust_obj, clust_assignment elements', call. = FALSE)
      
    }

    stopifnot(is_quosure(x$data))
    stopifnot(is.matrix(x$dist_mtx))
    stopifnot(x$dist.method %in% c(getDistMethods(), 'smc', 'sumofsquares'))
    stopifnot(x$clust_fun %in% c('hclust', 'kmeans', 'dbscan', 'som', 'pam', 'prediction'))
    stopifnot(is_tibble(x$clust_assignment))
    stopifnot(all(c('observation', 'clust_id') %in% names(x$clust_assignment)))
    
    structure(x, class = 'clust_analysis')

  }
  
  combi_analysis <- function(x) {
    
    ## creates a combi_analysis object
    ## given a result list
    
    if(!is.list(x)) stop('Only named lists can be converted to a clust_analysis objects', call. = FALSE)
    
    if(is.null(names(x))) stop('Only named lists can be converted to a clust_analysis objects', call. = FALSE)
    
    if(!all(c('clust_analyses', 'clust_assignment') %in% names(x))) {
      
      stop('The input list needs to provide clust_analyses and clust_assignment elements')
      
    }

    structure(x, class = 'combi_analysis')
    
  }
  
  red_analysis <- function(x) {
    
    ## creates a red_analysis object storing results of dimensinality reduction
    ## takes a list
    
    if(!is.list(x)) stop('Only named lists can be converted to a clust_analysis objects', call. = FALSE)
    
    if(is.null(names(x))) stop('Only named lists can be converted to a clust_analysis objects', call. = FALSE)
    
    if(!all(c('red_obj', 'red_fun', 'component_tbl', 'loadings') %in% names(x))) {
      
      stop('The input list needs to provide red_obj, red_fun, component_tbl and loadings elements')
      
    }
    
    structure(x, class = 'red_analysis')
    
  }
  
  dist <- function(x, ...) {
    
    UseMethod('dist', x)
    
  }
  
  dist.default <- function(x, ...) {
    
    stats::dist(x, ...)
    
  }
  
  extract <- function(x, ...) {
    
    UseMethod('extract', x)
    
  }
  
  components <- function(x, ...) {
    
    UseMethod('components', x)
    
  }
  
  var <- function(x, ...) {
    
    UseMethod('var', x)
    
  }
  
  var.default <- function(x, ...) {
    
    stats::var(x, ...)
    
  }
  
  nobs <- function(x, ...) {
    
    UseMethod('nobs', x)
    
  }
  
  nobs.default <- function(x, ...) {
    
    stats::nobs(x, ...)
    
  }
  
  ngroups <- function(x, ...) {
    
    UseMethod('ngroups', x)
    
  }
  
  cv <- function(x, ...) {
    
    UseMethod('cv', x)
    
  }
  
  impact <- function(x, ...) {
    
    UseMethod('impact', x)
    
  }
  
# OOP extractor and appearance methods ------
  
  extract.clust_analysis <- function(clust_analysis_object, 
                                     type = c('distance', 'assignment', 'clust_object', 'data', 'object')) {
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    type <- match.arg(type[1], 
                      c('distance', 'assignment', 'clust_object', 'data', 'object'))
    
    switch(type, 
           distance = dist(clust_analysis_object), 
           assignment = clust_analysis_object$clust_assignment, 
           clust_object = clust_analysis_object$clust_obj, 
           data = model.frame(clust_analysis_object), 
           object = clust_analysis_object$clust_obj)
    
  }
  
  extract.combi_analysis <- function(combi_analysis_object, 
                                     type = c('distance', 'assignment', 'clust_object', 'data', 'object')) {
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    if(type == 'assignment') {
      
      combi_analysis_object$clust_assignment
      
    } else {
      
      combi_analysis_object$clust_analyses %>% 
        map(extract, 
            type = type)
      
    }

  }
  
  extract.red_analysis <- function(red_analysis_object, 
                                   type = c('component_tbl', 'scores', 'loadings', 'data', 'sdev', 'object')) {
    
    stopifnot(class(red_analysis_object) == 'red_analysis')
    
    type <- match.arg(type[1], 
                      c('component_tbl', 'scores', 'loadings', 'data', 'sdev', 'object'))
    
    if(type != 'sdev') {
      
      switch(type, 
             component_tbl = red_analysis_object$component_tbl, 
             scores = red_analysis_object$component_tbl, 
             loadings = red_analysis_object$loadings, 
             data = model.frame(red_analysis_object), 
             object = red_analysis_object$red_obj)
      
    } else if(red_analysis_object$red_fun == 'pca') {
      
      tibble(component = 1:length(red_analysis_object$red_obj$sdev), 
             sdev = red_analysis_object$red_obj$sdev, 
             perc_sdev = red_analysis_object$red_obj$sdev/sum(red_analysis_object$red_obj$sdev)*100, 
             var = red_analysis_object$red_obj$sdev^2, 
             perc_var = red_analysis_object$red_obj$sdev^2/sum(red_analysis_object$red_obj$sdev^2)*100)
      
    } else {
      
      mds_dims <- extract(red_analysis_object, 'scores')
      
      dims <- stri_extract(names(mds_dims), regex = 'comp_\\d+')
      
      dims <- dims[!is.na(dims)]
  
      vars <- mds_dims[dims] %>% 
        map_dbl(var)
      
      sds <- mds_dims[dims] %>% 
        map_dbl(sd)
      
      tibble(component = 1:length(dims), 
             sdev = sds, 
             perc_sdev = sds/sum(sds)*100, 
             var = vars, 
             perc_var = vars/sum(vars)*100)
      
    }
    
  }
  
  model.frame.clust_analysis <- function(clust_analysis_object) {
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    eval_tidy(clust_analysis_object$data)
    
  }
  
  model.frame.red_analysis <- function(red_analysis_object) {
    
    stopifnot(class(red_analysis_object) == 'red_analysis')
    
    eval_tidy(red_analysis_object$data)
    
  }

  dist.clust_analysis <- function(clust_analysis_object) {
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    clust_analysis_object$dist_mtx
    
  }
  
  summary.clust_analysis <- function(clust_analysis_object) {
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    summary(clust_analysis_object$clust_obj)
    
  }
  
  print.clust_analysis <- function(clust_analysis_object) {
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    print(clust_analysis_object$clust_obj)
    print(clust_analysis_object$clust_assignment)
    
  }
  
  nobs.clust_analysis <- function(clust_analysis_object) {
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    extract(clust_analysis_object, 'data') %>% 
      get_data_dim

  }
  
  nobs.combi_analysis <- function(combi_analysis_object) {
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    extract(combi_analysis_object, 'data') %>% 
      get_data_dim
    
  }
  
  nobs.red_analysis <- function(red_analysis_object) {
    
    stopifnot(class(red_analysis_object) == 'red_analysis')
    
    extract(red_analysis_object, 'data') %>% 
      get_data_dim
    
  }
  
  ngroups.clust_analysis <- function(clust_analysis_object) {
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    extract(clust_analysis_object, 'assignment') %>% 
      count(clust_id)
    
  }
  
  ngroups.combi_analysis <- function(combi_analysis_object) {
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    extract(combi_analysis_object, 'assignment') %>% 
      count(clust_id)
    
  }
  
# OOP dimensionality reduction methods for the cluster objects ------

  components.clust_analysis <- function(clust_analysis_object, 
                                        kdim = NULL, 
                                        red_fun = c('pca', 'mds'), 
                                        with = c('distance', 'data'), ...) {
    
    ## performs k-dimensional MDS or PCA with the clust_analysis object
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    with <- match.arg(with[1], c('distance', 'data'))
    
    if(is.null(kdim)) {
      
      kdim <- length(unique(clust_analysis_object$clust_assignment$clust_id))
      
    }

    red_analysis_obj <- reduce_data(data = extract(clust_analysis_object, type = with), 
                                    distance_method = clust_analysis_object$dist_method, 
                                    kdim = kdim, 
                                    red_fun = red_fun, ...)
    
    red_analysis_obj$component_tbl <- red_analysis_obj$component_tbl %>% 
      left_join(extract(clust_analysis_object, 'assignment'), 
                by = 'observation')
    
    return(red_analysis_obj)
    
  }
  
  components.combi_analysis <- function(combi_analysis_object, 
                                        kdim = NULL, 
                                        red_fun = c('pca', 'mds'), 
                                        with = c('distance', 'data'), ...) {
    
    ## performs k-dimensional MDS or PCA with the observation-clustering combi_analysis object
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    with <- match.arg(with[1], c('distance', 'data'))
    
    red_analysis <- components(combi_analysis_object$clust_analyses$observation, 
                               kdim = kdim, 
                               red_fun = red_fun, 
                               with = with, ...)
    
    red_analysis$component_tbl <- red_analysis$component_tbl %>% 
      select(observation, 
             starts_with('comp')) %>% 
      left_join(extract(combi_analysis_object, 'assignment'), 
                by = 'observation')
    
    return(red_analysis)
    
  }
  
# OOP plotting methods ------
  
  plot.red_analysis <- function(red_analysis_object, 
                                type = c('component_tbl', 'scores', 'loadings', 'scree'), 
                                label_points = TRUE, 
                                cust_theme = theme_classic(), 
                                segment_color = 'steelblue', ...) {
    
    ## points the requested elements of the dimensionality reduction
    
    stopifnot(class(red_analysis_object) == 'red_analysis')
    
    if(red_analysis_object$red_fun == 'mds' & (type %in% c('loadings'))) {
      
      stop('Loadings are not implemented for the MDS dimansionality reduction algotithm', call. = FALSE)
      
    }
    
    type <- match.arg(type[1], 
                      c('component_tbl', 'scores', 'loadings', 'scree'))
    
    ## plot meta
    
    plot_n <- nobs(red_analysis_object)

    plot_tag <- paste0('\nObservations: n = ', 
                       plot_n$observations, 
                       '\nVariables: n = ', 
                       plot_n$variables)
    
    sdevs <- extract(red_analysis_object, 'sdev')
    
    if(red_analysis_object$red_fun == 'mds') {
      
      ax_labs <- map2(c('Dim 1', 'Dim 2'), 
                      signif(sdevs$perc_var[1:2], 3), 
                      ~paste0(.x, ', ', .y, '%'))
      
    } else {
      
      ax_labs <- map2(c('PC1', 'PC2'), 
                      signif(sdevs$perc_var[1:2], 3), 
                      ~paste0(.x, ', ', .y, '%'))
      
    }
    
    if(type %in% c('component_tbl', 'scores')) {
      
      plot_point(data = extract(red_analysis_object, 'scores'), 
                 x_var = 'comp_1', 
                 y_var = 'comp_2', 
                 fill_var = NULL, 
                 label_var = NULL, 
                 plot_title = switch(red_analysis_object$red_fun, 
                                     mds = 'MDS dimensions', 
                                     pca = 'PCA scores'), 
                 plot_tag = plot_tag, 
                 x_lab = ax_labs[[1]], 
                 y_lab = ax_labs[[2]], 
                 cust_theme = cust_theme, ...)
      
    } else if(type == 'loadings') {
      
      plot_point(data = extract(red_analysis_object, 'loadings'), 
                 x_var = 'comp_1', 
                 y_var = 'comp_2', 
                 fill_var = NULL, 
                 label_var = if(label_points) 'variable' else NULL, 
                 plot_title = 'PCA loadings', 
                 plot_tag = plot_tag, 
                 x_lab = ax_labs[[1]], 
                 y_lab = ax_labs[[2]], 
                 cust_theme = cust_theme, 
                 show_segments = TRUE, 
                 segment_color = segment_color, ...)
      
    } else {
      
      sdevs %>%
        mutate(line_group = 'gr1') %>% 
        ggplot(aes(x = component, 
                   y = perc_var, 
                   group = line_group)) + 
        geom_line(color = segment_color) + 
        cust_theme + 
        labs(title = switch(red_analysis_object$red_fun, 
                            mds = 'MDS variance', 
                            pca = 'PCA variance'), 
             y = '% total variance', 
             x = switch(red_analysis_object$red_fun, 
                        mds = 'Dimension', 
                        pca = 'Principal component'), 
             tag = plot_tag) +
        scale_x_continuous(breaks = sdevs$component)
      
    }
    
  }
  
  plot.clust_analysis <- function(clust_analysis_object, 
                                  type = c('diagnostic', 'components', 'heat_map', 'training'), 
                                  cust_theme = theme_classic(), 
                                  jitter_width = 0, 
                                  jitter_height = 0, 
                                  point_alpha = 1, ...) {
    
    ## plots diagnostic graphics
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    type <- match.arg(type[1], c('diagnostic', 'components', 'heat_map', 'training'))
    
    ## plot meta
    
    clust_algorithm <- clust_analysis_object$clust_fun %>% 
      switch(hclust = paste(',', clust_analysis_object$hc_method, 'method'), 
             kmeans = '', 
             pam = '', 
             som = '')
    
    clust_method <- clust_analysis_object$clust_fun %>% 
      switch(hclust = 'Hierarchical clustering', 
             kmeans = 'Kmeans clustering', 
             pam = 'PAM clustering', 
             som = 'SOM clustering', 
             prediction = 'prediction')
    
    plot_subtitle <- paste0(clust_method, 
                            ', ', 
                            clust_analysis_object$dist_method, 
                            ' distance', 
                            clust_algorithm)
    
    plot_n <- nobs(clust_analysis_object)
    
    plot_tag <- paste0('\nObservations: n = ', 
                       plot_n$observations, 
                       '\nVariables: n = ', 
                       plot_n$variables)
    
    ## plots
    
    if(type == 'diagnostic') {
      
      ## diagnostic plots
      
      plot_list <- list()
      
      if(clust_analysis_object$clust_fun %in% c('hclust', 'kmeans', 'pam')) {
        
        plot_title <- 'Optimal cluster number'
        
        k <- length(unique(clust_analysis_object$clust_assignment$clust_id))
        
        clust_args <- clust_analysis_object$clust_fun %>% 
          switch(hclust = list(FUNcluster = hcut, 
                               hc_method = clust_analysis_object$hc_method), 
                 kmeans = list(FUNcluster = hcut), 
                 pam = list(FUNcluster = pam))
        
        plot_list[c('wss', 'silhouette')] <- c('wss', 'silhouette') %>% 
          map(function(x) call2('plot_nbclust', 
                                dist_mtx = clust_analysis_object$dist_mtx,
                                method = x, 
                                k = k, 
                                !!!clust_args, 
                                plot_title = plot_title, 
                                plot_subtitle = plot_subtitle, 
                                plot_tag = plot_tag, 
                                cust_theme = cust_theme)) %>% 
          map(eval)
        
        if(clust_analysis_object$clust_fun == 'hclust') {
          
          plot_list$dendrogram <- plot_dendro(clust_str = clust_analysis_object$clust_obj, 
                                              k = k, 
                                              labels = T, 
                                              cust_theme = cust_theme, 
                                              plot_tag = plot_tag)
          
        }
        
      } else if(clust_analysis_object$clust_fun == 'dbscan') {
        
        plot_list$knn_dist <- plot_knn_distance(as.dist(clust_analysis_object$dist_mtx), 
                                                k = clust_analysis_object$minPts - 1, 
                                                eps = clust_analysis_object$eps, 
                                                plot_title = 'kNN distance plot', 
                                                plot_subtitle = paste('DB Scan algorithm, eps =', 
                                                                      clust_analysis_object$eps), 
                                                plot_tag = plot_tag, 
                                                cust_theme = cust_theme)
        
      } else if(clust_analysis_object$clust_fun == 'som') {
        
        plot_list <- plot_som(clust_analysis_object$clust_obj)
        
      } else {
        
        warning('No training plots available for cluster predictions', call. = FALSE)
        
        return(NULL)
        
      }
      
      return(plot_list)
      
    } else if(type == 'components') {
      
      ## MDS or PCA
      
      red_results <- components(clust_analysis_object, ...)
      
      stopifnot(all(c('comp_1', 'comp_2') %in% names(extract(red_results, 'scores'))))
      
      sdevs <- extract(red_results, 'sdev')
      
      if(red_results$red_fun == 'pca') {
        
        ax_labs <- map2(c('PC1', 'PC2'), 
                        signif(sdevs$perc_var[1:2], 3), 
                        ~paste0(.x, ', ', .y, '%'))
        
      } else {
        
        ax_labs <- map2(c('Dim 1', 'Dim 2'), 
                        signif(sdevs$perc_var[1:2], 3), 
                        ~paste0(.x, ', ', .y, '%'))
        
      }
      
      plot_point(data = extract(red_results, 'scores'), 
                 x_var = 'comp_1', 
                 y_var = 'comp_2', 
                 fill_var = 'clust_id', 
                 plot_title = switch(red_results$red_fun , 
                                     pca = 'PCA', 
                                     mds = 'MDS'), 
                 plot_subtitle = plot_subtitle, 
                 plot_tag = plot_tag, 
                 x_lab = ax_labs[[1]], 
                 y_lab = ax_labs[[2]], 
                 cust_theme = cust_theme, 
                 jitter_height = jitter_height, 
                 jitter_width = jitter_width, 
                 fill_lab = 'Cluster ID', 
                 point_alpha = point_alpha)
      
    } else if(type == 'heat_map') {
      
      ## heat map of the distances between the observations
      
      plotting_tbl <- clust_analysis_object$dist_mtx %>% 
        as.data.frame
      
      clust_vars <- colnames(plotting_tbl)
      
      plotting_tbl <- plotting_tbl %>% 
        rownames_to_column('observation') %>%
        left_join(clust_analysis_object$clust_assignment[c('observation', 'clust_id')], 
                  by = 'observation') %>% 
        gather(key = 'observation2', 
               value = 'distance', 
               all_of(clust_vars)) %>% 
        left_join(set_names(clust_analysis_object$clust_assignment[c('observation', 'clust_id')], 
                            c('observation2', 
                              'clust_id2')), 
                  by = 'observation2')
      
      heat_map <- plotting_tbl %>% 
        ggplot(aes(x = reorder(observation, distance), 
                   y = reorder(observation2, distance), 
                   fill = distance)) + 
        geom_tile() + 
        scale_fill_gradient2(low = 'firebrick', 
                             mid = 'white', 
                             high = 'steelblue', 
                             midpoint = mean(range(plotting_tbl$distance))) + 
        cust_theme + 
        theme(axis.title = element_blank(), 
              axis.text.x = element_text(angle = 90, 
                                         hjust = 1, 
                                         vjust = 0.5)) + 
        labs(title = 'Distances between observations', 
             subtitle = plot_subtitle, 
             tag = plot_tag,
             fill = 'Distance') + 
        facet_grid(clust_id2 ~ clust_id, 
                   scales = 'free', 
                   space = 'free', ...)
      
      return(heat_map)
      
    } else {
      
      if(clust_analysis_object$clust_fun != 'som') {
        
        warning('The training plots available only for the SOM cluster analyses', call. = FALSE)
        
        return(NULL)
        
      }
      
      som_training <- plot_train_som(kohonen_object = clust_analysis_object$clust_obj, 
                                     plot_title = 'SOM clustering: training', 
                                     plot_subtitle = plot_subtitle, 
                                     cust_theme = cust_theme)
      
      return(som_training)
      
    }
    
  }
  
  
  plot.combi_analysis <- function(combi_analysis_object, 
                                  type = c('diagnostic', 'components', 'heat_map', 'training'), 
                                  cust_theme = theme_classic(), 
                                  jitter_width = 0, 
                                  jitter_height = 0, 
                                  point_alpha = 1, ...) {
    
    ## plotting methods for the combi_analysis objects
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    type <- match.arg(type[1], c('diagnostic', 'components', 'heat_map', 'training'))
    
    ## plots for the step1 and step2 clustering analyses
    
    plot_list <- combi_analysis_object$clust_analyses %>% 
      map(plot, 
          type = type, 
          cust_theme = cust_theme, 
          jitter_height = jitter_height, 
          jitter_width = jitter_width, ...)
    
    ## summary component plots for the final clustering results
    
    if(type == 'components') {
      
      final_ass <- extract(combi_analysis_object, 'assignment')
      
      obs_red <- components(combi_analysis_object$clust_analyses$observation, ...)
      
      score_tbl <- extract(obs_red, 'scores') %>% 
        select(observation, 
               starts_with('comp')) %>% 
        left_join(final_ass, 
                  by = 'observation')
      
      sdevs <- extract(obs_red, 'sdev')
      
      if(obs_red$red_fun == 'pca') {
        
        ax_labs <- map2(c('PC1', 'PC2'), 
                        signif(sdevs$perc_var[1:2], 3), 
                        ~paste0(.x, ', ', .y, '%'))
        
      } else {
        
        ax_labs <- map2(c('Dim 1', 'Dim 2'), 
                        signif(sdevs$perc_var[1:2], 3), 
                        ~paste0(.x, ', ', .y, '%'))
        
      }
      
      plot_list$final <- plot_point(data = score_tbl, 
                                    x_var = 'comp_1', 
                                    y_var = 'comp_2', 
                                    fill_var = 'clust_id', 
                                    plot_title = switch(obs_red$red_fun , 
                                                        pca = 'PCA', 
                                                        mds = 'MDS'), 
                                    plot_tag = plot_list$observation$labels$tag, 
                                    x_lab = ax_labs[[1]], 
                                    y_lab = ax_labs[[2]], 
                                    cust_theme = cust_theme, 
                                    jitter_height = jitter_height, 
                                    jitter_width = jitter_width, 
                                    fill_lab = 'Cluster ID', 
                                    point_alpha = point_alpha)
      
    }
    
    return(plot_list)
    
  }
  
# OOP cluster quality control -----
  
  var.red_analysis <- function(red_analysis_object) {
    
    ## extracts component variances
    
    stopifnot(class(red_analysis_object) == 'red_analysis')
    
    extract(red_analysis_object, 'sdev')
    
  }
  
  var.clust_analysis <- function(clust_analysis_object) {
    
    ## calculates within and total sum of squares for any clustering object
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    get_sum_sq(dist_mtx = dist(clust_analysis_object), 
               assignment = extract(clust_analysis_object, 'assignment'))
    
  }
  
  var.combi_analysis <- function(combi_analysis_object) {
    
    ## calculates within and total sum of squares for any combi_analysis object
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    get_sum_sq(dist_mtx = extract(combi_analysis_object, 'distance')[[1]], 
               assignment = extract(combi_analysis_object, 'assignment'))
    
  }
  
# OOP semi-supervised clustering ------
  
  predict.clust_analysis <- function(clust_analysis_object, 
                                     newdata = NULL, 
                                     type = c('class', 'propagation'), ...) {
    
    ## assigns the observations from newdata to the clusters
    ## class: simple matching by observation names
    ## propagation: kNN-driven label propagation
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    type <- match.arg(type[1], 
                      c('class', 'propagation'))
    
    if(is.null(newdata)) {
      
      return(extract(clust_analysis_object, 'assignment'))
      
    }
    
    if(type == 'class') {
      
      train_assignment <- extract(clust_analysis_object, 'assignment') %>% 
        column_to_rownames('observation')
      
      if(nrow(newdata) != nrow(train_assignment)) {
        
        stop('The numbers of rows in new data and the table used for cluster development must be equal', call. = FALSE)
        
      }
      
      newdata <- as.data.frame(newdata)
      
      if(!is.null(rownames(newdata))) {
        
        test_assignment <- tibble(observation = rownames(newdata), 
                                  clust_id = train_assignment[rownames(newdata), 'clust_id'])
        
      } else {
        
        warning('Unnamed observations in new data')
        
        test_assignment <- train_assignment %>% 
          rownames_to_column('observation') %>%
          as_tibble
        
      }
      
      ## output
      
      model_frame <- enexpr(newdata)
      
      list(data = quo(model_frame), 
           dist_mtx = calculate_dist(newdata, method = clust_analysis_object$dist_method), 
           dist_method = clust_analysis_object$dist_method, 
           clust_fun = 'prediction', 
           clust_obj = NULL, 
           clust_assignment = test_assignment) %>% 
        clust_analysis
      
    } else {
      
      propagate(clust_analysis_object = clust_analysis_object, 
                newdata = newdata, 
                distance_method = clust_analysis_object$dist_method, ...)
      
    }
    
  }
  
  predict.combi_analysis <- function(combi_analysis_object, 
                                     newdata = NULL, 
                                     type = c('class', 'propagation'), ...) {
    
    ## assigns the observations from newdata to the clusters
    ## class: simple matching by observation names
    ## propagation: kNN-driven label propagation
    ## returns a cluster analysis object
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    type <- match.arg(type[1], 
                      c('class', 'propagation'))
    
    if(is.null(newdata)) {
      
      return(extract(combi_analysis_object, 'assignment'))
      
    }
    
    if(type == 'class') {
      
      train_assignment <- extract(combi_analysis_object, 'assignment') %>% 
        column_to_rownames('observation')
      
      if(nrow(newdata) != nrow(train_assignment)) {
        
        stop('The numbers of rows in new data and the table used for cluster development must be equal', call. = FALSE)
        
      }
      
      newdata <- as.data.frame(newdata)
      
      if(!is.null(rownames(newdata))) {
        
        test_assignment <- tibble(observation = rownames(newdata), 
                                  clust_id = train_assignment[rownames(newdata), 'clust_id'])
        
      } else {
        
        warning('Unnamed observations in new data')
        
        test_assignment <- train_assignment %>% 
          rownames_to_column('observation') %>%
          as_tibble
        
      }
      
      ## output
      
      model_frame <- enexpr(newdata)
      
      list(data = quo(model_frame), 
           dist_mtx = calculate_dist(newdata, method = combi_analysis_object$clust_analyses$observation$dist_method), 
           dist_method = combi_analysis_object$clust_analyses$observation$dist_method, 
           clust_fun = 'prediction', 
           clust_obj = NULL, 
           clust_assignment = test_assignment) %>% 
        clust_analysis
      
    } else {
      
      obs_propagation <- propagate(clust_analysis_object = combi_analysis_object$clust_analyses$observation, 
                                   newdata = newdata, 
                                   distance_method = combi_analysis_object$clust_analyses$observation$dist_method, ...)
      
      obs_propagation <- extract(obs_propagation, 'assignment') %>% 
        set_names(c('observation', 'node'))
      
      node_ass <- extract(combi_analysis_object, 'assignment') %>% 
        filter(!duplicated(node))

      test_assignment <- left_join(obs_propagation, 
                node_ass[c('node', 'clust_id')], 
                by = 'node')
      
      ## output
      
      model_frame <- enexpr(newdata)
      
      list(data = quo(model_frame), 
           dist_mtx = calculate_dist(newdata, method = combi_analysis_object$clust_analyses$observation$dist_method), 
           dist_method = combi_analysis_object$clust_analyses$observation$dist_method, 
           clust_fun = 'prediction', 
           clust_obj = NULL, 
           clust_assignment = test_assignment) %>% 
        clust_analysis
      
    }
    
  }

# OOP cluster stability -----
  
  cv.clust_analysis <- function(clust_analysis_object, 
                                nfolds = 5, 
                                nearest_n = 5, 
                                simple_vote = TRUE, 
                                kernel_fun = function(x) 1/x, 
                                seed = 1234, 
                                .parallel = FALSE) {
    
    ## cross validation of an existing clustering object
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    ## common parameters
    
    clust_data <- extract(clust_analysis_object, 'data')
    
    distance_method <- clust_analysis_object$dist_method
    
    if(clust_analysis_object$clust_fun == 'hclust') {
      
      cv_call <- call2('cv_cluster', 
                       data = clust_data, 
                       nfolds = nfolds, 
                       nearest_n = nearest_n, 
                       simple_vote = simple_vote, 
                       kernel_fun = kernel_fun,
                       distance_method = distance_method, 
                       clustering_fun = hcluster, 
                       seed = seed, 
                       .parallel = .parallel, 
                       k = nrow(ngroups(clust_analysis_object)), 
                       hc_method = clust_analysis_object$hc_method, 
                       !!!clust_analysis_object$dots)
      
    } else if(clust_analysis_object$clust_fun %in% c('kmeans', 'pam')) {
      
      cv_call <- call2('cv_cluster', 
                       data = clust_data, 
                       nfolds = nfolds, 
                       nearest_n = nearest_n, 
                       simple_vote = simple_vote, 
                       kernel_fun = kernel_fun,
                       distance_method = distance_method, 
                       clustering_fun = kcluster, 
                       clust_fun = clust_analysis_object$clust_fun, 
                       seed = seed, 
                       .parallel = .parallel, 
                       k = nrow(ngroups(clust_analysis_object)), 
                       !!!clust_analysis_object$dots)
      
    } else if(clust_analysis_object$clust_fun == 'som') {

        cv_call <- call2('cv_cluster', 
                         data = clust_data, 
                         nfolds = nfolds, 
                         nearest_n = nearest_n, 
                         simple_vote = simple_vote, 
                         kernel_fun = kernel_fun,
                         distance_method = distance_method, 
                         clustering_fun = som_cluster, 
                         seed = seed, 
                         .parallel = .parallel, 
                         xdim = clust_analysis_object$grid$xdim, 
                         ydim = clust_analysis_object$grid$ydim, 
                         topo = clust_analysis_object$grid$topo, 
                         neighbourhood.fct = as.character(clust_analysis_object$grid$neighbourhood.fct), 
                         toroidal = clust_analysis_object$grid$toroidal, 
                         !!!clust_analysis_object$dots)
      
    }
    
    eval(cv_call)
    
  }
  
  cv.combi_analysis <- function(combi_analysis_object, 
                                nfolds = 5, 
                                nearest_n = 5, 
                                simple_vote = TRUE, 
                                kernel_fun = function(x) 1/x, 
                                seed = 1234, 
                                .parallel = FALSE) {
    
    ## cross validation of an existing combi object
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    ## common parameters
    
    node_clust_fun <- switch(combi_analysis_object$clust_analyses$node$clust_fun, 
                             hclust = hcluster, 
                             kmeans = kcluster, 
                             pam = kcluster)
    
    cv_call <- call2('cv_cluster', 
                     data = extract(combi_analysis_object$clust_analyses$observation, 'data'), 
                     nfolds = nfolds, 
                     nearest_n = nearest_n, 
                     simple_vote = simple_vote, 
                     kernel_fun = kernel_fun, 
                     clustering_fun = combi_cluster, 
                     seed = seed, 
                     .parallel = .parallel, 
                     distance_som = combi_analysis_object$clust_analyses$observation$dist_method,
                     xdim = combi_analysis_object$clust_analyses$observation$grid$xdim, 
                     ydim = combi_analysis_object$clust_analyses$observation$grid$ydim, 
                     topo = combi_analysis_object$clust_analyses$observation$grid$topo, 
                     neighbourhood.fct = as.character(combi_analysis_object$clust_analyses$observation$grid$neighbourhood.fct), 
                     toroidal = combi_analysis_object$clust_analyses$observation$grid$toroidal, 
                     rlen = nrow(combi_analysis_object$clust_analyses$observation$clust_obj$changes), 
                     node_clust_fun = node_clust_fun, 
                     distance_nodes = combi_analysis_object$clust_analyses$node$dist_method, 
                     k = nrow(ngroups(combi_analysis_object$clust_analyses$node)), 
                     !!!combi_analysis_object$dots)
    
    eval(cv_call)
    
  }
  
# OOP variable importance -----
  
  impact.clust_analysis <- function(clust_analysis_object, 
                                    seed = 1234, 
                                    .parallel = FALSE) {
    
    ## gets impact of the clustering variables by noising
    
    stopifnot(class(clust_analysis_object) == 'clust_analysis')
    
    ## common parameters
    
    clust_data <- extract(clust_analysis_object, 'data')
    
    distance_method <- clust_analysis_object$dist_method
    
    if(clust_analysis_object$clust_fun == 'hclust') {
      
      imp_call <- call2('importance_cluster', 
                       data = clust_data, 
                       clustering_fun = hcluster, 
                       distance_method = distance_method, 
                       seed = seed, 
                       .parallel = .parallel, 
                       k = nrow(ngroups(clust_analysis_object)), 
                       hc_method = clust_analysis_object$hc_method, 
                       !!!clust_analysis_object$dots)
      
    } else if(clust_analysis_object$clust_fun %in% c('kmeans', 'pam')) {
      
      imp_call <- call2('importance_cluster', 
                        data = clust_data, 
                        clustering_fun = kcluster, 
                        clust_fun = clust_analysis_object$clust_fun,
                        distance_method = distance_method, 
                        seed = seed, 
                        .parallel = .parallel, 
                        k = nrow(ngroups(clust_analysis_object)), 
                        !!!clust_analysis_object$dots)
      

    } else if(clust_analysis_object$clust_fun == 'som') {
      
      imp_call <- call2('importance_cluster', 
                        data = clust_data, 
                        clustering_fun = hcluster, 
                        clust_fun = som_cluster,
                        distance_method = distance_method, 
                        seed = seed, 
                        .parallel = .parallel, 
                        xdim = clust_analysis_object$grid$xdim, 
                        ydim = clust_analysis_object$grid$ydim, 
                        topo = clust_analysis_object$grid$topo, 
                        neighbourhood.fct = as.character(clust_analysis_object$grid$neighbourhood.fct), 
                        toroidal = clust_analysis_object$grid$toroidal,
                        !!!clust_analysis_object$dots)

    }
    
    eval(imp_call)
    
  }
  
  impact.combi_analysis <- function(combi_analysis_object, 
                                    seed = 1234, 
                                    .parallel = FALSE) {
    
    ## gets impact of the clustering variables by noising
    
    stopifnot(class(combi_analysis_object) == 'combi_analysis')
    
    ## common paramaters
    
    node_clust_fun <- switch(combi_analysis_object$clust_analyses$node$clust_fun, 
                             hclust = hcluster, 
                             kmeans = kcluster, 
                             pam = kcluster)
    
    imp_call <- call2('importance_cluster', 
                      data = extract(combi_analysis_object$clust_analyses$observation, 'data'), 
                      clustering_fun = combi_cluster, 
                      seed = seed, 
                      .parallel = .parallel, 
                      distance_som = combi_analysis_object$clust_analyses$observation$dist_method,
                      xdim = combi_analysis_object$clust_analyses$observation$grid$xdim, 
                      ydim = combi_analysis_object$clust_analyses$observation$grid$ydim, 
                      topo = combi_analysis_object$clust_analyses$observation$grid$topo, 
                      neighbourhood.fct = as.character(combi_analysis_object$clust_analyses$observation$grid$neighbourhood.fct), 
                      toroidal = combi_analysis_object$clust_analyses$observation$grid$toroidal, 
                      rlen = nrow(combi_analysis_object$clust_analyses$observation$clust_obj$changes), 
                      node_clust_fun = node_clust_fun, 
                      distance_nodes = combi_analysis_object$clust_analyses$node$dist_method, 
                      k = nrow(ngroups(combi_analysis_object$clust_analyses$node)), 
                      !!!combi_analysis_object$dots)
    
    eval(imp_call)
    
    
  }
  
# hierarchical, kmeans, pam, dbscan and som clustering -----
  
  get_clust_tendency <- function(data, n, seed = 1234, ...) {
    
    ## checks clustering tendency with Hopkins statistic
    ## a wrapper around the respective function provided by factoextra
    
    tend <- factoextra::get_clust_tendency(data = data, 
                                           n = n, 
                                           seed = seed, ...)

    tend$p_value <- 1 - pbeta(tend$hopkins_stat, 
                              shape1 = n, 
                              shape2 = n)
    
    return(tend)
        
  }
  
  hcluster <- function(data, 
                       distance_method = 'euclidean', 
                       k = 2, 
                       hc_method = 'ward.D2', 
                       seed = 1234, ...) {
    
    ## calculates distance between the data in the data frame with named rows using philentropy package
    ## the simple matching distance id provided by scrime
    ## and subjects the distance matrix to clustering using the hierarchical clustering algorithm
    ## provides plot of mean squared error vs cluster number and silhouette stat vs cluster number 
    ## to facilitate the decision on the cluster number
    
    set.seed(seed = seed)
    
    ## distance calculation
    
    dist_mtx <- calculate_dist(data = data, 
                               method = distance_method)
    
    ## hierarchical clustering and cluster assignment table
    
    hclust_str <- hclust(dist_mtx %>% 
                           as.dist, 
                         method = hc_method, ...)
    
    hclust_ass <- hclust_str %>% 
      cutree(k = k)
    
    hclust_ass <- tibble(observation = names(hclust_ass), 
                         clust_id = factor(unname(hclust_ass)))
    
    ## output
    
    model_frame <- enexpr(data)
    
    list(data = quo(model_frame), 
         dist_mtx = dist_mtx, 
         dist_method = distance_method, 
         clust_fun = 'hclust', 
         clust_obj = hclust_str, 
         clust_assignment = hclust_ass, 
         hc_method = hc_method, 
         dots = list2(...)) %>% 
      clust_analysis
    
  }
  
  kcluster <- function(data, 
                       distance_method = 'euclidean', 
                       clust_fun = 'kmeans',
                       k = 2, 
                       seed = 1234, ...) {
    
    ## calculates distance between the data in the data frame with named rows using philentropy package
    ## the simple matching distance is provided by scrime
    ## and subjects the distance matrix to clustering using the kmeans or pam clustering algorithm
    ## provides plot of mean squared error vs cluster number and silhouette stat vs cluster number 
    ## to facilitate the decision on the cluster number
    
    set.seed(seed = seed)
    
    ## distance calculation
    
    dist_mtx <- calculate_dist(data = data, 
                               method = distance_method) %>% 
      as.dist
    
    ## kmeans/pam clustering and cluster assignment table
    
    arguments <- list(dist_mtx, k, ...)
    
    kclust_str <- call2(clust_fun, !!!arguments) %>% 
      eval

    kclust_ass <- tibble(observation = names(kclust_str$cluster), 
                         clust_id = factor(unname(kclust_str$cluster)))

    ## output
    
    model_frame <- enexpr(data)

    list(data = quo(model_frame), 
         dist_mtx = as.matrix(dist_mtx), 
         dist_method = distance_method, 
         clust_fun = clust_fun, 
         clust_obj = kclust_str, 
         clust_assignment = kclust_ass, 
         dots = list2(...)) %>% 
      clust_analysis
    
  }
  
  dbscan_cluster <- function(data, 
                             distance_method = 'euclidean', 
                             eps, 
                             minPts = 5, ...) {
    
    ## calculates distance between the data in the data frame with named rows using philentropy package, 
    ## the simple matching distance is provided by scrime
    ## and subjects the distance matrix to clustering using the dbscan clustering algorithm
    ## to facilitate the decision on the optimal eps value, a plot of kNN distance (where k = minPts - 1)
    ## versus sample is provided
    
    ## distance calculation
    
    dist_mtx <- calculate_dist(data = data, 
                               method = distance_method)

    ## dbscan clustering and cluster assignment table
    
    kclust_str <- dbscan(dist_mtx %>% 
                           as.dist, 
                         eps = eps, 
                         minPts = minPts, ...)
    
    kclust_ass <- tibble(observation = rownames(dist_mtx), 
                         clust_id = factor(kclust_str$cluster))
    
    
    ## output
    
    model_frame <- enexpr(data)
    
    list(data = quo(model_frame), 
         dist_mtx = dist_mtx, 
         dist_method = distance_method, 
         clust_fun = 'dbscan', 
         clust_obj = kclust_str, 
         clust_assignment = kclust_ass,
         eps = eps, 
         minPts = minPts, 
         dots = list2(...)) %>% 
      clust_analysis
    
  }
  
  som_cluster <- function(data, 
                          distance_method = 'euclidean',
                          xdim = floor(sqrt(5*sqrt(nrow(data)))), 
                          ydim = floor(sqrt(5*sqrt(nrow(data)))), 
                          topo = 'hexagonal', 
                          neighbourhood.fct = 'gaussian', 
                          toroidal = FALSE, 
                          seed = 1234, ...) {
    
    ## performs SOM clustering of the data
    ## warning: still limited number of distance kernels implemented
    
    ## grid
    
    som_grid <- somgrid(xdim = xdim, 
                        ydim = ydim, 
                        topo = topo, 
                        neighbourhood.fct = neighbourhood.fct, 
                        toroidal = toroidal)

    ## fitting

    set.seed(seed = seed)
    
    kohonen_obj <- som(as.matrix(data), 
                       grid = som_grid, 
                       dist.fcts = distance_method, ...)
    
    node_ass <- tibble(observation = rownames(kohonen_obj$data[[1]]), 
                       clust_id = factor(kohonen_obj$unit.classif), 
                       node = factor(kohonen_obj$unit.classif), 
                       neuro_dist = kohonen_obj$distances)
    
    ## output
    
    model_frame <- enexpr(data)
    
    list(data = quo(model_frame), 
         dist_mtx = calculate_dist(data = data, 
                                   method = distance_method), 
         dist_method = distance_method, 
         clust_fun = 'som', 
         clust_obj = kohonen_obj, 
         clust_assignment = node_ass,
         grid = som_grid, 
         dots = list2(...)) %>% 
      clust_analysis

  }
  
# combined SOM clustering -----
  
  combi_cluster <- function(data, 
                            distance_som = 'euclidean',
                            xdim = floor(sqrt(5*sqrt(nrow(data)))), 
                            ydim = floor(sqrt(5*sqrt(nrow(data)))), 
                            topo = 'hexagonal', 
                            neighbourhood.fct = 'gaussian', 
                            toroidal = FALSE, 
                            rlen = NULL, 
                            node_clust_fun = hcluster, 
                            distance_nodes = 'euclidean', 
                            k = 3, 
                            seed = 1234, ...) {
    
    ## a combi wrapper for two-step clustering of the study data:
    ## Step 1: som with the provided grid object 
    ## Step 2: hierarchical/k-mean like clustering with the k-branches/centers
    
    ## SOM clustring of the observations
    
    clust_analyses <- list()
    
    som_clust <- som_cluster(data = data, 
                             distance_method = distance_som, 
                             xdim = xdim, 
                             ydim = ydim, 
                             topo = topo, 
                             neighbourhood.fct = neighbourhood.fct, 
                             toroidal = toroidal, 
                             rlen = rlen, 
                             seed = seed)
    
    ## second-step clustering
    
    node_clust <- node_clust_fun(data = som_clust$clust_obj$codes[[1]], 
                                 distance_method = distance_nodes, 
                                 k = k, 
                                 seed = seed, ...)
    
    ## assignment to the nodes and node clusters
    
    som_ass <- extract(som_clust, 'assignment') %>% 
      select(observation, node)
    
    node_ass <- extract(node_clust, 'assignment') %>% 
      select(observation, clust_id) %>% 
      mutate(observation = stri_replace(observation, fixed = 'V', replacement = '')) %>% 
      set_names(c('node', 'clust_id'))

    combi_ass <- left_join(som_ass, 
                           node_ass, 
                           by = 'node') 
    
    ## output
    
    list(clust_analyses = list(observation = som_clust, 
                               node = node_clust), 
         clust_assignment = combi_ass, 
         dots = list2(...)) %>% 
      combi_analysis
    
  }
  
# dimensionality reduction -----
  
  reduce_data <- function(data, 
                          distance_method = 'euclidean', 
                          kdim = 2, 
                          red_fun = c('pca', 'mds'), ...) {
    
    ## performs dimensionality reduction of the provided data using MDS or PCA
    ## the distance method is only valid for MDS
    
    red_fun <- match.arg(red_fun[1], c('pca', 'mds'))
    
    if(red_fun == 'mds') {
      
      red_obj <- NULL
      
      dist_mtx <- calculate_dist(data = data, 
                                 method = distance_method)
      
      component_tbl <- dist_mtx %>% 
        as.dist %>% 
        cmdscale(k = kdim, ...)
      
      loadings <- NULL
      
    } else {
      
      red_obj <- PCAproj(x = data, 
                         k = kdim, ...)
      
      component_tbl <- red_obj$scores
      
      rownames(component_tbl) <- rownames(data)
      
      loadings <- red_obj$loadings %>% 
        unclass %>%  
        as.data.frame %>% 
        rownames_to_column('variable') %>% 
        as_tibble
      
      loadings <- set_names(loadings, 
                            c('variable', 
                              paste('comp', 1:(ncol(loadings) - 1), sep = '_')))
      
    }

    component_tbl <- component_tbl %>% 
      as.data.frame %>% 
      mutate(observation = rownames(data))  %>% 
      as_tibble
    
    component_tbl <- set_names(component_tbl, 
                               c(paste('comp', 1:(ncol(component_tbl) - 1), sep = '_')), 
                               'observation')
    
    ## output
    
    model_frame <- enexpr(data)
      
    list(data = quo(model_frame), 
         red_obj = red_obj, 
         red_fun = red_fun, 
         component_tbl = component_tbl, 
         loadings = loadings) %>% 
      red_analysis
    
  }
  
# semi-supervised learning -----
  
  propagate <- function(clust_analysis_object, 
                        newdata = NULL, 
                        distance_method = clust_analysis_object$dist_method, 
                        k = 5, 
                        simple_vote = TRUE, 
                        kernel_fun = function(x) 1/x, 
                        detailed = FALSE) {
    
    ## predicts cluster assignment by k-NN driven label propagation
    ## simple_vote: classical kNN prediction, ties are resolved by returning the first element otherwise, distance weighting
    ## detailed: kNN calculation results (distance, ids and labels are returned as well)
    ## value: a clust_analysis object
    
    if(class(clust_analysis_object) != 'clust_analysis') {
      
      stop('A valid clust_analysis object required')
      
    }
    
    if(is.null(newdata)) {
      
      return(extract(clust_analysis_object, 'assignment'))
      
    }

    ## extracting the training data set, constrained to the cases with the cluster assignment (non-NA)
    
    train_set <- extract(clust_analysis_object, type = 'data') %>% 
      as.matrix
    
    if(ncol(train_set) != ncol(newdata)) {
      
      stop('The numbers of columns in new data and the table used for cluster development must be equal', call. = FALSE)
      
    }
    
    if(is.null(rownames(train_set))) {
      
      rownames(train_set) <- paste0('train_', 1:nrow(train_set))
      
    }
    
    train_assign <- extract(clust_analysis_object, type = 'assignment') %>% 
      mutate(.rowname = rownames(train_set)) %>% 
      filter(!is.na(clust_id))
    
    label_vec <- set_names(train_assign$clust_id, 
                           train_assign$.rowname)
    
    train_set <- train_set[train_assign$.rowname, ]
    
    ## constructing the mixed test-train matrix and calculating a distance object
    
    if(is.null(rownames(newdata))) {
      
      rownames(newdata) <- paste0('test_', 1:nrow(train_set))
      
    }
    
    mix_set <- rbind(train_set, as.matrix(newdata))
    
    mix_diss <- calculate_dist(data = mix_set, 
                               method = distance_method)
    
    ## kNN calulation and label assingment
    
    knn_dists <- kNN(as.dist(mix_diss), k = nrow(mix_diss) - 1)
    
    knn_test <- list(dist = knn_dists$dist[rownames(newdata), ], 
                     id = knn_dists$id[rownames(newdata), ])
    
    knn_test$annot_id <- matrix(rownames(mix_diss)[knn_test$id], ncol = nrow(mix_diss) - 1)
    rownames(knn_test$annot_id) <- rownames(knn_test$id)
    
    knn_test$labels <- matrix(label_vec[knn_test$annot_id], ncol = nrow(mix_diss) - 1)
    rownames(knn_test$labels) <- rownames(knn_test$id)
    
    ## voting, constrained to the k-nearest neighbors from the train dataset
    
    if(simple_vote) {
      
      clust_assignment <- rownames(knn_test$labels) %>% 
        map(~knn_test$labels[.x, ]) %>% 
        map(~.x[!is.na(.x)][1:k]) %>% 
        map(vote_simple) %>% 
        unlist
      
    } else {
      
      rows <- rownames(knn_test$labels) %>% 
        map(~knn_test$labels[.x, ])
      
      dists <- rownames(knn_test$labels) %>% 
        map(~knn_test$dist[.x, ])
      
      non_na <- rows %>% 
        map(~!is.na(.x))
      
      rows <- map2(rows, non_na, ~.x[.y][1:k])
      dists <- map2(dists, non_na, ~.x[.y][1:k])
      
      clust_assignment <- list(vector = rows, 
                               dist_vec = dists) %>% 
        pmap(vote_kernel, 
             kernel_fun = kernel_fun) %>% 
        unlist

    }

    model_frame <- enexpr(newdata)
    
    clust_analysis_out <- list(data = quo(model_frame), 
                               dist_mtx = as.matrix(mix_diss)[rownames(newdata), rownames(newdata)], 
                               dist_method = distance_method, 
                               clust_fun = 'prediction', 
                               clust_obj = NULL, 
                               clust_assignment = tibble(observation = rownames(knn_test$labels), 
                                                         clust_id = factor(clust_assignment))) %>% 
      clust_analysis
    
    if(!detailed) {
      
      return(clust_analysis_out)
      
    } else {
      
      return(list(mix_data = mix_set, 
                  mix_diss = mix_diss, 
                  knn = knn_test[c('dist', 'id', 'labels')], 
                  clust_analysis_object = clust_analysis_out))
      
    }
    
  }
  
# cross-validation of the clustering alghorithm stability -----
  
  cv_cluster <- function(data,
                         nfolds = 5, 
                         nearest_n = 5, 
                         simple_vote = TRUE, 
                         kernel_fun = function(x) 1/x, 
                         clustering_fun = kcluster, 
                         seed = 1234, 
                         .parallel = FALSE, ...) {
    
    ## tests stability of the clustering algorithm by cross validation
    
    start_time <- Sys.time()
    message(paste('CV for', nfolds, 'folds'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    set.seed(seed = seed)
    
    ## data split generation
    
    fold_ids <- createFolds(1:nrow(data), k = nfolds)

    fold_set <- fold_ids %>% 
      map(~list(train = data[-.x, ], 
                test = data[.x, ])) %>% 
      transpose
    
    ## defining the global classifier
    
    glob_classif <- clustering_fun(data = data, ...)
    
    glob_assign <- extract(glob_classif, 'assignment')[c('observation', 'clust_id')] %>%
      set_names(c('observation', 'global_clust'))
    
    ## creating the training set classifiers
    
    if(.parallel) {
      
      plan('multisession')

      train_classif <- fold_set$train %>% 
        future_map(function(x) clustering_fun(data = x, ...), 
                   .options = furrr_options(seed = TRUE, 
                                            packages = c('tidyverse', 
                                                         'rlang', 
                                                         'cluster', 
                                                         'kohonen', 
                                                         'dbscan', 
                                                         'Rcpp', 
                                                         'somKernels'), 
                                            globals = c('hcluster', 
                                                        'kcluster', 
                                                        'som_cluster', 
                                                        'dbscan_cluster', 
                                                        'combi_cluster', 
                                                        'clust_analysis', 
                                                        'combi_analysis', 
                                                        'calculate_dist', 
                                                        'getDistMethods', 
                                                        'source_kernels', 
                                                        'extract', 
                                                        'extract.clust_analysis', 
                                                        'stri_replace', 
                                                        'sm')))
      
      plan('sequential')
      
    } else {
      
      train_classif <- fold_set$train %>% 
        map(function(x) clustering_fun(data = x, ...))
      
    }
    
    ## obtaining the predictions and comparing with the global classifier
    
    if(class(glob_classif) == 'clust_analysis') {
      
      pred_input <- list(clust_analysis_object = train_classif, 
                         newdata = fold_set$test)
      
    } else {
      
      pred_input <- list(combi_analysis_object = train_classif, 
                         newdata = fold_set$test)
      
    }
    
    test_preds <- pred_input %>% 
      pmap(safely(predict), 
           type = 'propagation', 
           k = nearest_n, 
           simple_vote = simple_vote, 
           kernel_fun = kernel_fun) %>% 
      map(~.x$result) %>% 
      compact %>% 
      map(extract, 'assignment')
    
    test_preds <- test_preds %>% 
      map(select, observation, clust_id) %>% 
      map(set_names, c('observation', 'fold_clust')) %>% 
      map(left_join, glob_assign, by = 'observation') %>% 
      map(mutate, correct = as.character(global_clust) == as.character(fold_clust))

    test_stats <- test_preds %>% 
      map(summarise, 
          corr_rate = mean(as.numeric(correct), na.rm = TRUE)) %>% 
      map2_dfr(., names(.), 
               ~mutate(.x, err_rate = 1 - corr_rate, fold = .y, )) %>% 
      select(fold, corr_rate, err_rate)
    
    bca_err <- bca(test_stats$err_rate)
    
    test_summary <- tibble(mean_error = mean(test_stats$err_rate, na.rm = TRUE), 
                           lower_ci = bca_err[1], 
                           upper_ci = bca_err[2])
    
    test_preds <- test_preds %>% 
      map2_dfr(., names(.), ~mutate(.x, fold = .y))
    
    ## output
    
    list(clust_analysis_object = glob_classif, 
         predictions = test_preds, 
         fold_stats = test_stats, 
         summary = test_summary)
    
  }
  
# variable importance determined by noising -------
  
  importance_cluster <- function(data, 
                                 clustering_fun = kcluster, 
                                 seed = 1234, 
                                 .parallel = FALSE, ...) {
    
    ## determines importance of the clustering variables by noising
    ## as proposed by Leo Breiman for random forests
    
    ## noised set and the unnoised data
    
    noised_set <- names(data) %>% 
      map(function(x) mutate(data, !!ensym(x) := sample(.data[[x]], size = nrow(data), replace = TRUE))) %>% 
      set_names(names(data))
    
    noised_set <- c(list(data = data), noised_set)
    
    ## generating the clustering objects, calculating the variances
    
    if(.parallel) {
      
      plan('multisession')
      
      var_lst <- noised_set %>% 
        future_map(function(x) clustering_fun(data = x, ...), 
                   .options = furrr_options(seed = TRUE, 
                                            packages = c('tidyverse', 
                                                         'rlang', 
                                                         'cluster', 
                                                         'kohonen', 
                                                         'dbscan', 
                                                         'Rcpp', 
                                                         'somKernels'), 
                                            globals = c('hcluster', 
                                                        'kcluster', 
                                                        'som_cluster', 
                                                        'dbscan_cluster', 
                                                        'combi_cluster', 
                                                        'clust_analysis', 
                                                        'combi_analysis', 
                                                        'calculate_dist', 
                                                        'getDistMethods', 
                                                        'source_kernels', 
                                                        'extract', 
                                                        'extract.clust_analysis', 
                                                        'stri_replace', 
                                                        'sm')))
      
      plan('sequential')
      
    } else {
      
      var_lst <- noised_set %>% 
        map(function(x) clustering_fun(data = x, ...))
      
    }
    
    var_lst <- var_lst %>% 
      map(var) %>% 
      map(~.x[c('total_wss', 'total_ss', 'between_ss', 'frac_var')]) %>% 
      map_dfr(as_tibble) %>% 
      mutate(variable = c('data', names(data)))
    
    var_summary <- var_lst[-1, ] %>% 
      mutate(frac_diff = var_lst$frac_var[1] - frac_var)
    
    list(variances = var_lst, 
         summary = var_summary[c('variable', 'frac_diff')])
    
  }
  
# Varia -----
  
  plot_clust_hm <- function(sample_clust_object, 
                            feature_clust_object = NULL, 
                            plot_title = NULL, 
                            plot_subtitle = NULL, 
                            x_lab = 'sample', 
                            cust_theme = theme_classic(), 
                            discrete_fill = FALSE) {
    
    ## plots clustering features as a heat map
    
    if(!class(sample_clust_object) %in% c('combi_analysis', 'clust_analysis')) {
      
      stop('The sample cluster object has to be of clust_analysis or combi_analysis class', call. = FALSE)
      
    }
    
    if(!is.null(feature_clust_object)) {
      
      if(!class(feature_clust_object) %in% c('combi_analysis', 'clust_analysis')) {
        
        stop('The sample cluster object has to be of clust_analysis or combi_analysis class', call. = FALSE)
        
      }
      
    }

    ## assignment list
    
    if('node' %in% names(sample_clust_object$clust_assignment)) {
      
      sample_ass_names <- c('sample', 'sample_node', 'sample_clust')
      
    } else {
      
      sample_ass_names <- c('sample', 'sample_clust')
      
    }
  
    if(!is.null(feature_clust_object)) {
      
      feature_ass_names <- switch(class(feature_clust_object), 
                                  clust_analysis = c('feature', 'feature_clust'), 
                                  combi_analysis = c('feature', 'feature_node', 'feature_clust'))
      
      cmm_assignment <- list(sample_clust_object, 
                             feature_clust_object) %>% 
        map(extract, type = 'assignment') %>% 
        map2(list(sample_ass_names, 
                  feature_ass_names), 
             set_names) %>% 
        set_names(c('sample', 'feature'))
      
    } else {
      
      cmm_assignment <- extract(sample_clust_object, 'assignment') %>% 
        set_names(sample_ass_names) %>% 
        list(sample = .)
      
    }
    
    ## data in a long format
    
    data <- switch(class(sample_clust_object), 
                   clust_analysis = extract(sample_clust_object, 'data'), 
                   combi_analysis = extract(sample_clust_object, 'data')[[1]]) %>% 
      as_tibble 
      
    features <- names(data)
    
    data <- data %>% 
      mutate(sample = cmm_assignment$sample$sample) %>% 
      gather(key = 'feature', 
             value = 'value', 
             all_of(features))
    
    ## joining the data table with the cluster assignment information
    
    data <- left_join(data, 
                      cmm_assignment$sample, 
                      by = 'sample')

    if(!is.null(feature_clust_object)) {
      
      if(all(!features %in% cmm_assignment$feature$feature)) {
        
        stop('Unable to retrieve the cluster assignment information for the clustering variables from the feature_clust_object.')
        
      }
      
      data <- left_join(data, 
                        cmm_assignment$feature, 
                        by = 'feature')
      
    }

    ## plot tag with the sample n numbers
    
    n_tag <- ngroups(sample_clust_object)
    
    n_tag <- map2_chr(n_tag$clust_id, 
                      n_tag$n, 
                      ~paste0(.x, ': n = ', .y)) %>% 
      paste(collapse = ', ')
    
    ## plotting
    
    if(class(sample_clust_object) == 'clust_analysis') {
      
      base_plot <- data %>% 
        ggplot(aes(x = reorder(sample, as.numeric(factor(value))), 
                   y = reorder(feature, as.numeric(factor(value))), 
                   fill = if(discrete_fill) factor(value) else value))
      
    } else {
      
      base_plot <- data %>% 
        ggplot(aes(x = reorder(sample, as.numeric(sample_node)), 
                   y = reorder(feature, as.numeric(factor(value))), 
                   fill = if(discrete_fill) factor(value) else value))
      
    }
    
    if(!is.null(feature_clust_object)) {
      
      base_plot <- base_plot + 
        facet_grid(feature_clust ~ sample_clust, 
                   scales = 'free', 
                   space = 'free')
      
    } else {
      
      base_plot <- base_plot + 
        facet_grid(. ~ sample_clust, 
                   scales = 'free', 
                   space = 'free')
      
    }
    
    base_plot +
      geom_tile() + 
      cust_theme + 
      theme(axis.title.y = element_blank(), 
            panel.background = element_blank(), 
            axis.line = element_blank(), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank()) + 
      labs(title = plot_title,
           subtitle = plot_subtitle, 
           tag = n_tag, 
           x = x_lab)

  }
  
# END ----