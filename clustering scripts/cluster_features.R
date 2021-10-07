# This script compares the clustering features between the clinical clusters

  insert_head()
  
# data container -----
  
  clust_features <- list()
  
# globals -----
  
  ## variable colors and labels
  
  clust_features$var_labs <- c('Age at diagnosis', 
                               'Female sex', 
                               'SMWD', 
                               'WHO func. class', 
                               '5-year mortality') %>% 
    set_names(study_clust$variables)
  
  clust_features$var_colors <- c('gray60', 
                                 'indianred', 
                                 'chartreuse4', 
                                 'steelblue3', 
                                 'coral4') %>% 
    set_names(study_clust$variables)

# serial analysis -----
  
  insert_msg('Serial analysis')
  
  ## analysis objects
  
  clust_features$analyses <- study_clust$analysis_tbl %>% 
    map(function(cohort) study_clust$variables %>% 
          map(analyze_feature, 
              inp_tbl = cohort, 
              split_var = 'clust_id') %>% 
          set_names(study_clust$variables)) %>% 
    set_names(names(study_clust$analysis_tbl))
  
  ## summaries
  
  clust_features$summaries <- clust_features$analyses %>% 
    map(function(cohort) cohort %>% 
          map_dfr(extract_test_summary) %>% 
          filter(test == 'anova') %>% 
          mutate(p_adj = p.adjust(p_value, 'BH')))
  
  ## counts
  
  clust_features$stats <- clust_features$analyses %>% 
    map(function(cohort) cohort %>% 
          map_dfr(extract_counts))
  
# Displaying the feature quantities as radial plots -----
  
  insert_msg('Displaying the feature quantities as radial plots')
  
  clust_features$plots <- list(clust_stats = clust_features$stats, 
                               plot_title = globals$center_labs[names(clust_features$stats)]) %>% 
    pmap(plot_radial) %>% 
    map(function(x) x + 
          scale_fill_manual(values = clust_features$var_colors, 
                            labels = clust_features$var_labs, 
                            name = '') + 
          facet_grid(. ~ split_var, 
                     labeller = as_labeller(study_clust$clust_labs)))
  
# END -----
  
  insert_tail()