# Principal component analysis of the study features following min/max normalization

  insert_head()
  
# container list -----
  
  pca <- list()
  
# globals: analysis tables -----
  
  insert_msg('Globals setup')
  
  pca$variables <- pah_study$mod_variables$variable
  
  pca$analysis_tbl <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(select, ID, all_of(pca$variables)) %>% 
    map(~map_dfc(.x, function(var) if(is.factor(var)) as.numeric(var) else var)) %>% 
    map(~map_dfc(.x, function(var) if(is.numeric(var)) scale(var, center = median(var), scale = TRUE)[, 1] else var)) %>% 
    map(column_to_rownames, 'ID')
  
# PCA dimension number ------
  
  insert_msg('Dimension number for PCA')
  
  pca$dim_no_pca <- pca$analysis_tbl %>% 
    map(reduce_data, 
        kdim = 8, 
        red_fun = 'pca')
  
  ## scree plots: 4 dimensions explain 90% data
  
  pca$dim_no_scree <- pca$dim_no_pca %>% 
    map(plot, 
        type = 'scree', 
        cust_theme = globals$common_theme)
  
# 4-dimensional PCA ------
  
  insert_msg('4D PCA')
  
  pca$pca_obj <- pca$analysis_tbl %>% 
    map(reduce_data, 
        kdim = 4, 
        red_fun = 'pca')
  
  pca$pca_obj$IBK_0$loadings <- pca$pca_obj$IBK_0$loadings %>% 
    mutate(variable = translate_vars(variable))
  
  pca$pca_obj$LZ_0$loadings <- pca$pca_obj$LZ_0$loadings %>% 
    mutate(variable = translate_vars(variable))
  
  ## Plots of PCA scores and loadings
  
  pca$pca_score_plots <- list(red_analysis_object = pca$pca_obj,  
                              point_color = globals$center_colors[c('IBK_0', 'LZ_0')]) %>% 
    pmap(plot, 
        type = 'scores', 
        cust_theme = globals$common_theme, 
        plot_subtitle = 'PCA scores') %>% 
    map2(., globals$center_labs[c('IBK_0', 'LZ_0')], 
         ~.x + labs(title = .y))
  
  pca$pca_loadings_plots <- list(red_analysis_object = pca$pca_obj,  
                                 point_color = globals$center_colors[c('IBK_0', 'LZ_0')]) %>% 
    pmap(plot, 
         type = 'loadings', 
         cust_theme = globals$common_theme, 
         plot_subtitle = 'PCA loadings') %>% 
    map2(., globals$center_labs[c('IBK_0', 'LZ_0')], 
         ~.x + labs(title = .y))
  
# Plotting the eigenvector lengths as a point plot -----
  
  insert_msg('Eigenvector length')
  
  ## length table
  
  pca$eigen_length <- pca$pca_obj %>% 
    map(extract, 'loadings') %>% 
    map(mutate, 
        vec_len = comp_1^2 +comp_2^2 + comp_3^2 + comp_4^2, 
        vec_len = sqrt(vec_len))
  
  ## plot
  
  pca$eigen_plot <- list(data = pca$eigen_length, 
                         color = globals$center_colors[c('IBK_0', 'LZ_0')], 
                         title = globals$center_labs[c('IBK_0', 'LZ_0')]) %>% 
    pmap(function(data, color, title) ggplot(data, 
                                             aes(x = vec_len, 
                                                 y = reorder(variable, vec_len))) + 
           geom_point(shape = 16, 
                      size = 2, 
                      color = color) + 
           globals$common_theme + 
           theme(axis.title.y = element_blank()) + 
           labs(title = title, 
                subtitle = 'Factor influence', 
                x = 'Loadings vector length'))
  
# Assesing the general clustering tendency of the raw data and PCA scores ----
  
  insert_msg('Clustering tendency')
  
  pca$clust_tendency_data <- pca$analysis_tbl %>% 
    map(get_clust_tendency, 
        n = 50, 
        seed = 1234)
  
  pca$clust_tendency_pca <- pca$pca_obj %>% 
    map(extract, 'scores') %>% 
    map(column_to_rownames, 'observation') %>% 
    map(get_clust_tendency, 
        n = 50, 
        seed = 1234)
  
  
# END ------
  
  insert_msg()


  