# This script tries to identify subsets of participants differing in the following clinical features:
# (1) age as a readout of frailty
# (2) sex
# (3) SMWD as a readout of mobility
# (4) WHO functional class as readout of handicap
# (5) 5-year mortality as a readout of therapy success

  c('./tools/sys_tools.R', 
    './tools/clust_tools.R', 
    './tools/counting_tools.R', 
    './tools/project_tools.R') %>% 
    walk(source)

  insert_head()
  
# data container -----
  
  study_clust <- list()
  
# globals: analysis table -----
  
  insert_msg('Globals setup')
  
  ## variables
  
  study_clust$variables <- c('age_fc', 
                             'sex', 
                             'SMWD', 
                             'WHOFc',
                             'death_acute')
  
  ## analysis table: males coded as 0, min-max normalization
  
  study_clust$analysis_tbl <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(select, 
        all_of(study_clust$variables)) %>% 
    map(mutate, 
        sex = as.numeric(sex) - 1) %>% 
    map(function(cohort) cohort %>% 
          map_dfc(min_max)) %>% 
    map2(.,
         map(pah_study[c('IBK_0', 'LZ_0')], 
             ~.x[['ID']]), 
         set_rownames)
  
  ## cluster labs and colors
  
  study_clust$clust_labs <- c('cluster_1' = 'Sub 1', 
                              'cluster_2' = 'Sub 2', 
                              'cluster_3' = 'Sub 3', 
                              'noise' = 'UA')
  
  study_clust$clust_colors <- c('cluster_1' = 'steelblue', 
                                'cluster_2' = 'cornsilk4', 
                                'cluster_3' = 'indianred3', 
                                'noise' = 'lightgoldenrod3')
  
# DB scan clustering ----
  
  insert_msg('Clustering with DBSCAN algorithm')
  
  ## serial analysis
  
  set.seed(1234)
  
  study_clust$clust_analysis <- list(inp_tbl = study_clust$analysis_tbl,
                                     eps = c(0.57, 0.8)) %>% 
    pmap(dbscan_data, 
         distance_method = 'euclidean', 
         minPts = 5)
  
  ## assignment 
  
  study_clust$clust_analysis$IBK_0$clust_assignment <- study_clust$clust_analysis$IBK_0$clust_assignment %>% 
    mutate(clust_id = car::recode(clust_id, "1 = 'cluster_1'; 2 = 'cluster_2'; 3 = 'cluster_3'; 4 = 'noise'") %>% 
             factor)
  
  study_clust$clust_analysis$LZ_0$clust_assignment <- study_clust$clust_analysis$LZ_0$clust_assignment %>% 
    mutate(clust_id = car::recode(clust_id, "2 = 'cluster_1'; 1 = 'cluster_2'; 3 = 'cluster_3'; 4 = 'noise'") %>% 
             factor)
  
  ## n numbers
  
  study_clust$n_numbers <- study_clust$clust_analysis %>% 
    map(~.x[['clust_assignment']]) %>% 
    map(count, 
        clust_id) %>% 
    map(mutate, 
        clust_lab = study_clust$clust_labs[clust_id], 
        plot_tag = paste0(clust_lab, 
                          ': n = ', 
                          n))
  
  ## plot tags
  
  study_clust$plot_tags <- study_clust$n_numbers %>% 
    map(~.x[['plot_tag']]) %>% 
    map(paste, collapse = ', ') %>% 
    map(function(x) paste0('\n', x))
  
# Displaying the cluster assignment in PCA plots ----  
  
  insert_msg('Cluster PCA')
  
  ## base plots
  
  study_clust$pca_plots <- list(cluster_analysis = study_clust$clust_analysis, 
                                plot_title = paste('Clustering of participants', 
                                                   globals$center_labs[c('IBK_0', 'LZ_0')], 
                                                   sep = ': ')) %>% 
    pmap(plot_clust_mds, 
         k_dim = 2, 
         red_fun = 'pca', 
         cluster_labs = study_clust$clust_labs, 
         cluster_colors = study_clust$clust_colors, 
         cust_theme = globals$common_theme, 
         plot_subtitle = 'DB scan algorithm') %>% 
    transpose
  
  ## extracting the PCA variances
  
  study_clust$pca_plots$var_perc <- study_clust$pca_plots$red_obj %>% 
    map(~.x$sdev) %>%
    map(~.x^2) %>% 
    map(function(x) x/sum(x))
  
  ## updating the plots
  
  study_clust$pca_plots <- list(plot = study_clust$pca_plots$plot, 
                                x_lab = map(study_clust$pca_plots$var_perc, ~.x[[1]]), 
                                y_lab = map(study_clust$pca_plots$var_perc, ~.x[[2]]), 
                                plot_tag = study_clust$plot_tags) %>% 
    pmap(function(plot, x_lab, y_lab, plot_tag) plot + 
           labs(x = paste0('PC1, ', 
                          signif(x_lab * 100, 2), 
                          '%'), 
                y = paste0('PC2, ', 
                           signif(y_lab * 100, 2), 
                           '%'), 
                tag = plot_tag))
  
# Generating tables with the cluster assignment, clustering variables and PAH risk scores -----
  
  insert_msg('Generating tables with the clustering info, variables and risk scores')
  
  ## cluster assignment
  
  study_clust$analysis_tbl <- study_clust$analysis_tbl %>%
    map(rownames_to_column,
        'ID') %>% 
    map2(., 
         map(study_clust$clust_analysis, ~.x$clust_assignment) %>% 
           map(set_names, 
               c('ID', 'clust_id')), 
         left_join, by = 'ID')
  
  ## scores and comparators
  
  study_clust$analysis_tbl <- map2(study_clust$analysis_tbl, 
                                   multi_modeling$score_tbl, 
                                   cbind) %>% 
    map2(., 
         map(pah_study[c('IBK_0', 'LZ_0')], 
             select,
             all_of(pah_study$comparators$variable)), 
         cbind) %>% 
    map(as_tibble)
  
# Downstream analyses -----
  
  insert_msg('Downstream analyses')
  
  c('./clustering scripts/cluster_features.R', 
    './clustering scripts/cluster_overall_surv.R', 
    './clustering scripts/cluster_scores.R', 
    './clustering scripts/gender_survival.R') %>% 
    walk(source)
  
# END -----
  
  insert_tail()