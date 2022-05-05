# This script generates main and supplementary figures -----

  insert_head()

# data containers ----

  paper_figures <- list()
  suppl_figures <- list()
  
# Figure 1: STROBE ----
  
  insert_msg('Figure 1: STROBE')
  
  paper_figures$strobe <- plot_grid(ggdraw() + 
                                      draw_image('./aux files/strobe.png')) %>% 
    as_figure(label = 'figure_1_strobe', 
                     w = 180, 
                     h = 140)
  
# Figure 2: development of the elastic net model -----
  
  insert_msg('Figure 2: elastic net model')
  
  paper_figures$elanet$top_panel <- plot_grid(multi_plots$lasso_hr_plot + 
                                                theme(legend.position = 'none'), 
                                              ggdraw(), 
                                              ncol = 2, 
                                              rel_widths = c(0.7, 0.3))
  
  paper_figures$elanet$bottom_panel <- multi_plots$km_plots %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = .,
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(multi_plots$km_plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1))
  
  paper_figures$elanet <- plot_grid(paper_figures$elanet$top_panel, 
                                    paper_figures$elanet$bottom_panel, 
                                    nrow = 2, 
                                    labels = LETTERS, 
                                    label_size = 10, 
                                    rel_heights = c(0.54, 0.46)) %>% 
    as_figure(label = 'figure_2_elanet', 
              w = 180, 
              h = 200)
  
# Figure 3: model calibration ------   
  
  insert_msg('Figure 3: Calibration of the elastic net model')
  
  paper_figures$elanet_calibration <- multi_plots$km_quart_plots %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(multi_plots$km_quart_plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1)) %>% 
    as_figure(label = 'figure_3_elanet_tertiles', 
              w = 180,
              h = 90)
  
# Figure 4: clustering of the participants -----
  
  insert_msg('Figure 4: participant clustering')

  paper_figures$part_clust$top_panel <- clust$pca_plots %>% 
    map(~.x + 
          theme(legend.position = 'none', 
                plot.subtitle = element_blank(), 
                plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(clust$pca_plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1))
  
  paper_figures$part_clust$bottom_panel <- clust$vio_panels %>% 
    map(~.x + 
          theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(clust$vio_panels[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1))
  
  paper_figures$part_clust <- plot_grid(paper_figures$part_clust$top_panel, 
                                        paper_figures$part_clust$bottom_panel, 
                                        nrow = 2, 
                                        rel_heights = c(0.35, 0.65), 
                                        labels = LETTERS, 
                                        label_size = 10) %>% 
    as_figure(label = 'figure_4_cluster_development', 
              w = 180, 
              h = 180)
    
# Figure 5: differences in risk between the clusters -----
  
  insert_msg('Figure 5: differences in the the risk scales and survival between the clusters')
  
  paper_figures$cluster_risk$top_panel <- map2(cl_chara$plots[c('IBK_0.FRENCH3p', 
                                                                'IBK_0.FRENCH4p', 
                                                                'IBK_0.SPAHR', 
                                                                'IBK_0.Compera', 
                                                                'IBK_0.mRASP')], 
                                               cl_chara$plots[c('LZ_0.FRENCH3p', 
                                                                'LZ_0.FRENCH4p', 
                                                                'LZ_0.SPAHR', 
                                                                'LZ_0.Compera', 
                                                                'LZ_0.mRASP')], 
                                               set_common_y) %>% 
    map(~map(.x, ~.x + theme(plot.tag = element_blank(), 
                             legend.position = 'none', 
                             axis.text.x = element_blank(), 
                             axis.title.x = element_blank()))) %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 2, 
                   align = 'hv')) %>% 
    map(~.x + theme(plot.margin = ggplot2::margin(r = 3, l = 3, t = 2, unit = 'mm'))) %>% 
    c(list(get_legend(cl_chara$plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              labels = c('A', 'B', 'C', 'D', 'E'), 
              label_size = 10)
  
  paper_figures$cluster_risk$bottom_panel <- cl_surv$plots %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(cl_surv$plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1))
  
  paper_figures$cluster_risk <- plot_grid(paper_figures$cluster_risk$top_panel, 
                                          paper_figures$cluster_risk$bottom_panel, 
                                          nrow = 2, 
                                          rel_heights = c(0.6, 0.4), 
                                          labels = c('', 'F'), 
                                          label_size = 10) %>% 
    as_figure(label = 'figure_5_cluster_risk', 
              w = 180, 
              h = 210)
  
# Figure 6: analysis result summary -----
  
  insert_msg('Analysis result summary')
  
  paper_figures$summary <- ggdraw() + 
    draw_image('./aux files/summary_figure.png')
  
  paper_figures$summary <- paper_figures$summary %>% 
    as_figure(label = 'figure_6_summary', 
              w = 180, 
              h = 100)
  
# Figure S1: univariable Cox modeling ------
  
  insert_msg('Unvariate Cox modeling results')
  
  suppl_figures$uni_cox <- uni_cox$forest_plot + 
    theme(legend.position = 'bottom')
  
  suppl_figures$uni_cox <- suppl_figures$uni_cox %>% 
    as_figure(label = 'figure_s1_uni_cox', 
              w = 180, 
              h = 210)

# Figure S2: Cluster QC -----
  
  insert_msg('QC of the clustering object')
  
  suppl_figures$cluster_qc$top_panel <- plot_grid(clust_dev$test_plot + 
                                                    theme(legend.position = 'right'), 
                                                  ggdraw(), 
                                                  ncol = 2, 
                                                  rel_widths = c(0.8, 0.2), 
                                                  labels = c('A', ''), 
                                                  label_size = 10)
  
  suppl_figures$cluster_qc$bottom_panel <- plot_grid(clust$diagn_plots$wss + 
                                                       theme(plot.tag = element_blank()), 
                                                     clust$impact_plot + 
                                                       theme(plot.subtitle = element_blank()), 
                                                     ncol = 2, 
                                                     align = 'hv', 
                                                     labels = c('B', 'C'), 
                                                     label_size = 10)
    
  suppl_figures$cluster_qc <- plot_grid(suppl_figures$cluster_qc$top_panel, 
                                        suppl_figures$cluster_qc$bottom_panel, 
                                        nrow = 2, 
                                        rel_heights = c(0.55, 0.45)) %>% 
    as_figure(label = 'figure_s2_cluster_qc', 
              w = 180, 
              h = 190)
  
# Figure S3: Study features, differences in the clusters ------
  
  insert_msg('Differences in the study features between the clusters')
  
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
    as_figure(label = 'figure_s3_cluster_differences', 
              w = 180, 
              h = 120)
  
# saving the figures -----
  
  insert_msg('Saving the figures')
  
  paper_figures %>% 
    walk(save_figure, 
         path = './paper/figures', 
         device = cairo_pdf)
  
  paper_figures %>% 
    walk(save_figure, 
         format = 'png', 
         path = './paper/figures', 
         dpi = 600)
  
  suppl_figures %>% 
    walk(save_figure, 
         path = './paper/supplementary figures', 
         device = cairo_pdf)
  
# END ----
  
  insert_tail() 