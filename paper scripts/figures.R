# This script generates main and supplementary figures -----

  insert_head()

# data containers ----

  paper_figures <- list()
  suppl_figures <- list()
  rev_figures <- list()
  
# Figure 1: STROBE ----
  
  insert_msg('Figure 1: STROBE')
  
  paper_figures$strobe <- plot_grid(ggdraw() + 
                                      draw_image('./aux files/strobe.png')) %>% 
    as_figure(label = 'figure_1_strobe', 
              ref_name = 'strobe', 
              caption = paste('Flow diagram of the study', 
                              'analysis inclusion process.'), 
              w = 180, 
              h = 180 * 2166/3906)
  
# Figure 2: development of the elastic net model -----
  
  insert_msg('Figure 2: elastic net model')
  
  paper_figures$elanet$top_panel <- 
    plot_grid(multi_plots$coef_plot + 
                theme(legend.position = 'none'), 
              ggdraw(), 
              ncol = 2, 
              rel_widths = c(0.7, 0.3))
  
  paper_figures$elanet$bottom_panel <- multi_plots$km_quart_plots %>% 
    map(~.x + 
          labs(color = 'Elastic Net score', 
               subtitle = stri_replace(.x$labels$tag, 
                                       fixed = '\n', 
                                       replacement = '')) + 
          theme(legend.title = element_blank(), 
                legend.position = 'bottom', 
                plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  paper_figures$elanet <- plot_grid(paper_figures$elanet$top_panel, 
                                    paper_figures$elanet$bottom_panel, 
                                    nrow = 2, 
                                    labels = LETTERS, 
                                    label_size = 10, 
                                    rel_heights = c(0.54, 0.46)) %>% 
    as_figure(label = 'figure_2_elanet', 
              ref_name = 'elanet',
              caption = paste('Multi-parameter modeling', 
                              'of PAH survival with Elastic Net', 
                              'Cox regression.'), 
              w = 180, 
              h = 200)
  
# Figure 3: risk assessment tool comparison ------   
  
  insert_msg('Figure 3: Risk tool performance comparison')
  
  paper_figures$risk_tools <- surv_tools$fit_plots
  
  ## the RF ensemble is presented only in a figure for the Reviewer 2
  
  for(i in names(paper_figures$risk_tools)) {
    
    paper_figures$risk_tools[[i]]$data <- 
      paper_figures$risk_tools[[i]]$data %>% 
      filter(!stri_detect(variable, fixed = 'RF'))
    
  }
  
  paper_figures$risk_tools <- paper_figures$risk_tools %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr', 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_3_risk_tool_performance', 
              ref_name = 'risk_tools', 
              caption = paste('Performance of PAH risk assessment tools.'), 
              w = 180, 
              h = 90)

# Figure 4: clustering of the participants -----
  
  insert_msg('Figure 4: participant clustering')
  
  ## upper panel: UMAP layouts

  paper_figures$part_clust$top_panel <- clust$umap_plots %>% 
    map(~.x + 
          labs(x = stri_replace(.x$labels$x, 
                                fixed = 'Dim', 
                                replacement = 'UMAP'), 
               y = stri_replace(.x$labels$y, 
                                fixed = 'Dim', 
                                replacement = 'UMAP')) + 
          theme(legend.position = 'bottom',
                legend.title = element_blank(), 
                plot.subtitle = element_blank(), 
                plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv')
  
  ## bottom panel: violin plots, indicating significances in the Y axes
  
  paper_figures$part_clust$labels <- cl_chara$test %>% 
    map(filter, variable %in% clust$variables) %>% 
    map(mutate, 
        variable = exchange(variable, 
                           dict = globals$var_labs), 
        var_lab = paste(variable, significance, sep = '\n')) %>% 
    map(~set_names(.x$var_lab, .x$variable))
  
  paper_figures$part_clust$bottom_panel <- clust$vio_panels %>% 
    map2(., paper_figures$part_clust$labels, 
         ~.x + scale_y_discrete(labels = .y)) %>% 
    map(~.x + 
          labs(subtitle = .x$labels$tag %>% 
                 stri_replace(fixed = '\n', 
                              replacement = '')) + 
          theme(legend.position = 'none', 
                plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(clust$vio_panels[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1))
  
  paper_figures$part_clust <- 
    plot_grid(paper_figures$part_clust$top_panel, 
              paper_figures$part_clust$bottom_panel, 
              nrow = 2, 
              rel_heights = c(0.42, 0.57), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_4_cluster_development', 
              ref_name = 'part_clust', 
              caption = 'Clustering of the study participants.', 
              w = 180, 
              h = 190)
    
# Figure 5: differences in risk between the clusters -----
  
  insert_msg('Figure 5: risk scales and survival in the clusters')

  ## top panel: risk scales
  
  paper_figures$cluster_risk$top_panel <- cl_chara$plots %>% 
    map(~.x[c('FRENCH3p', 'FRENCH4p', 
              'SPAHR', 'Compera', 
              'mRASP')]) %>% 
    transpose
  
  paper_figures$cluster_risk$top_panel <- 
    map2(paper_figures$cluster_risk$top_panel, 
         list(c('steelblue4', 'steelblue2', 'coral2', 'coral4'), 
              c('steelblue4', 'steelblue2', 'coral2', 'coral4', 'gray60'), 
              c('steelblue', 'coral2', 'coral4'), 
              c('steelblue', 'coral2', 'coral4'), 
              c('steelblue', 'coral2', 'coral4')), 
         function(lst, colors) lst %>% 
           map(~.x + scale_fill_manual(values = colors))) %>% 
    unlist(recursive = FALSE)
  
  paper_figures$cluster_risk$top_panel <- 
    paper_figures$cluster_risk$top_panel %>% 
    map(~.x + 
          theme(legend.position = 'none', 
                axis.title.x = element_blank(), 
                plot.subtitle = element_blank())) %>% 
    c(list(get_legend(paper_figures$cluster_risk$top_panel$FRENCH4p.IBK_0 + 
                        labs(fill = 'FPHR, # risk factors'))), 
      list(get_legend(paper_figures$cluster_risk$top_panel$SPAHR.IBK_0 + 
                        labs(fill = 'SPAHR/COMPERA/mRASP\nrisk strata')))) %>% 
    plot_grid(plotlist = ., 
              ncol = 4, 
              align = 'hv', 
              axis = 'tblr', 
              labels = c('A', '', 'B', '', 
                         'C', '', 'D', '', 
                         'E', ''), 
              label_size = 10)

  ## bottom panel: Kaplan-Meier plots
  
  paper_figures$cluster_risk$bottom_panel <- cl_surv$plots %>% 
    map(~.x + 
          theme(legend.position = 'bottom', 
                legend.title = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv')
  
  paper_figures$cluster_risk <- 
    plot_grid(paper_figures$cluster_risk$top_panel, 
              paper_figures$cluster_risk$bottom_panel, 
              nrow = 2, 
              rel_heights = c(0.65, 0.35), 
              labels = c('', 'F'), 
              label_size = 10) %>% 
    as_figure(label = 'figure_5_cluster_risk', 
              ref_name = 'cluster_risk', 
              caption = paste('Risk assessment and survival differences', 
                              'in the participant clusters.'), 
              w = 180, 
              h = 220)
  
# Figure 6: analysis result summary -----
  
  insert_msg('Figure 6: Analysis result summary')
  
  paper_figures$summary <- ggdraw() + 
    draw_image('./aux files/summary_figure.png')
  
  paper_figures$summary <- paper_figures$summary %>% 
    as_figure(label = 'figure_6_summary', 
              ref_name = 'summary', 
              caption = 'Summary of the analysis results.', 
              w = 180, 
              h = 100)
  
# Figure S1: univariable Cox modeling ------
  
  insert_msg('Figure S1: Unvariate Cox modeling results')
  
  suppl_figures$uni_cox <- uni_cox$forest_plot + 
    theme(legend.position = 'bottom', 
          plot.tag = element_blank())
  
  suppl_figures$uni_cox <- suppl_figures$uni_cox %>% 
    as_figure(label = 'figure_s1_uni_cox', 
              ref_name = 'uni_cox', 
              caption = 'Univariable Cox proportional hazard modeling.', 
              w = 180, 
              h = 220)

# Figure S2: correlation of the Elastic Net score and risk assessment tools -----
  
  insert_msg('Figure S2: correlation of risk assessment tools')
  
  suppl_figures$tool_correlation <- corr_tools$plots %>% 
    map(~.x + 
          theme(legend.position = 'none', 
                axis.text.x = element_text(angle = 90, 
                                           hjust = 1, 
                                           vjust = 0.5))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    plot_grid(get_legend(corr_tools$plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1)) %>% 
    as_figure(label = 'figure_s2_risk_tool_correlation', 
              ref_name = 'tool_correlation', 
              caption = paste('Correlation of the Elastic Net score and', 
                              'established risk assessment scores.'), 
              w = 180, 
              h = 120)
  
# Figure S3: development of the LASSO ensemble model -----
  
  insert_msg('Figure S3: model ensemble')
  
  suppl_figures$ensemble$top_panel <- 
    plot_grid(lasso_tools$coef_plot + 
                theme(legend.position = 'none'), 
              ggdraw(), 
              ncol = 2, 
              rel_widths = c(0.7, 0.3))
  
  suppl_figures$ensemble$bottom_panel <- surv_tools$calibration_plots %>% 
    map(~.x$ensemble_score) %>% 
    map(~.x + 
          labs(color = 'Elastic Net score', 
               subtitle = stri_replace(.x$labels$tag, 
                                       fixed = '\n', 
                                       replacement = '')) + 
          theme(legend.title = element_blank(), 
                legend.position = 'bottom', 
                plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')
  
  suppl_figures$ensemble <- plot_grid(suppl_figures$ensemble$top_panel, 
                                      suppl_figures$ensemble$bottom_panel, 
                                      nrow = 2, 
                                      labels = LETTERS, 
                                      label_size = 10, 
                                      rel_heights = c(0.54, 0.46)) %>% 
    as_figure(label = 'figure_s3_ensemble', 
              ref_name = 'ensemble',
              caption = paste('Development of a LASSO ensemble model', 
                              'including the Elastic Net score and established', 
                              'PAH risk assessment tools.'), 
              w = 180, 
              h = 200)
  
# Figure S4: Cluster QC -----
  
  insert_msg('Figure S4: QC of the clustering object')
  
  ## top panel: searching for the optimal algorithm
  
  suppl_figures$cluster_qc$top_panel <- 
    plot_grid(clust_dev$test_plot + 
                theme(legend.position = 'right'), 
              ggdraw(), 
              ncol = 2, 
              rel_widths = c(0.9, 0.1))
  
  ## bottom panel: optimal cluster number
  
  suppl_figures$cluster_qc$bottom_panel <- 
    c(clust$diagn_plots[c("wss", "silhouette")], 
         clust_dev$k_plots["k_cv"]) %>% 
    map(~.x + theme(legend.position = 'none', 
                    plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr')

  suppl_figures$cluster_qc <- 
    plot_grid(suppl_figures$cluster_qc$top_panel, 
              suppl_figures$cluster_qc$bottom_panel, 
              nrow = 2, 
              rel_heights = c(0.4, 0.6), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s4_cluster_qc', 
              ref_name = 'cluster_qc', 
              caption = 'Development of the PAH clusters.', 
              w = 180, 
              h = 210)
  
# Figure S5: clusterng factor importance -------
  
  insert_msg('Figure S6: clustering factor importance')

  suppl_figures$importance <- clust$impact_plot %>% 
    as_figure(label = 'figure_s5_importance', 
              ref_name = 'importance', 
              caption = paste('Permutation importance of the variables', 
                              'used for development of the PAH clusters.'), 
              w = 110, 
              h = 110)
    
# Figure S6: Study features, differences in the clusters ------
  
  insert_msg('Figure S6: Differences in the study features between the clusters')
  
  suppl_figures$clust_diff <- cl_chara$summary_plots %>% 
    map(~.x + 
          scale_x_continuous(limits = c(0, 0.8)) +  
          theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              align = 'hv', 
              ncol = 2) %>% 
    plot_grid(get_legend(cl_chara$summary_plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1)) %>%
    as_figure(label = 'figure_s6_cluster_differences', 
              ref_name = 'clust_diff', 
              caption = paste('Differences in study variables between', 
                              'the PAH clusters.'), 
              w = 180, 
              h = 120)
  
# Figure R1: clustering with the pooled dataset --------
  
  insert_msg('Figure R1: clustering with the pooled dataset')
  
  rev_figures$pool_clust <- 
    plot_grid(clust_rev$test_plot + 
                theme(legend.position = 'right')) %>% 
    as_figure(label = 'figure_r1_pooled_data_clustering', 
              ref_name = 'pool_clust', 
              caption = paste('Clustering analysis with the', 
                              'pooled IBK and LZ/W data set.'), 
              w = 180, 
              h = 110)
  
# Figure R2: Random Forest ensemble --------
  
  insert_msg('Figure R2: Random Forest ensemble')
  
  rev_figures$rf_ensemble <- surv_tools$fit_plots %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(plot_grid(rf_tools$importance$plot, 
                        ncol = 2), 
              ., 
              nrow = 2, 
              rel_heights = c(0.7, 1), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_r2_random_forest_ensemble', 
              ref_name = 'rf_ensemble', 
              caption = paste('Development of a Random Forest ensemble of', 
                              'PAH risk assessment tools and its performance', 
                              'at predicting overall survival.'), 
              w = 180, 
              h = 150)
  
# saving the figures -----
  
  insert_msg('Saving the figures')
  
  paper_figures %>% 
    walk(pickle, 
         path = './paper/figures', 
         device = cairo_pdf)
  
  #paper_figures %>% 
   # walk(pickle, 
    #     format = 'png', 
      #   path = './paper/figures', 
      #   dpi = 600)
  
  suppl_figures %>% 
    walk(pickle, 
         path = './paper/supplementary figures', 
         device = cairo_pdf)
  
  rev_figures %>% 
    walk(pickle, 
         path = './paper/reviewer figures', 
         device = cairo_pdf)
  
# END ----
  
  rm(i)
  
  insert_tail() 