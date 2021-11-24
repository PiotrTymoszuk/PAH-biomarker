# This script generates main and supplementary figures -----

  insert_head()

# data containers ----

  paper_figures <- list()
  suppl_figures <- list()
  
# Figure 1: STROBE ----
  
  insert_msg('Figure 1: STROBE')
  
  paper_figures$strobe <- plot_grid(ggdraw() + 
                                      draw_image('./aux files/strobe.png')) %>% 
    as_figure_object(figure_label = 'figure_1_strobe', 
                     w = 180, 
                     h = 150)
  
# Figure 2: results of univariate Cox modeling ------
  
  insert_msg('Figure 2: univariate Cox modeling')
  
  paper_figures$uni_cox <- uni_modeling$forest_plot + 
    labs(subtitle = 'Univariable Cox proporional hazard modeling') + 
    theme(legend.position = 'bottom')
  
  paper_figures$uni_cox <- paper_figures$uni_cox %>% 
    as_figure_object(figure_label = 'figure_2_univariate_cox', 
                     w = 180, 
                     h = 180)
  
# Figure 3: performance of the developed risk signatures at predicting the overall survival ----
  
  insert_msg('Figure 3: predictive performance of the developed risk signatures, overall survival')
  
  paper_figures$signatures_os <- plot_grid(ggdraw() + 
                                             draw_image('./aux files/signature_os.png') + 
                                             theme(plot.margin = globals$common_margin), 
                                           multi_plots$c_plots$testing_c + 
                                             labs(subtitle = 'Cox proportional hazard modeling'), 
                                           nrow = 2,
                                           rel_heights = c(0.18, 0.82), 
                                           labels = LETTERS, 
                                           label_size = 10) %>% 
    as_figure_object(figure_label = 'figure_3_signatures_os', 
                     w = 180, 
                     h = 210)
  
# Figure 4: 5-year survival -----
  
  insert_msg('Figure 4: predictive performance of the developed risk signatures, 5-year survival')
  
  paper_figures$signatures_five <- plot_grid(ggdraw() + 
                                               draw_image('./aux files/signature_five.png') + 
                                               theme(plot.margin = globals$common_margin), 
                                             multi_plots$five_surv_logis + 
                                               labs(subtitle = 'Logistic modeling') + 
                                               theme(legend.position = 'bottom'), 
                                             nrow = 2, 
                                             rel_heights = c(0.15, 0.85), 
                                             labels = LETTERS, 
                                             label_size = 10) %>% 
    as_figure_object(figure_label = 'figure_4_signatures_five', 
                     w = 180, 
                     h = 230)
  
# Figure 5: Gender and age-specific survival differences -----
  
  insert_msg('Figure 5: gender/age interplay and survival')
  
  ## the top panel
  
  paper_figures$gender_surv$upper_panel <- gender$os$plots[c('IBK_young', 
                                                             'LZ_young', 
                                                             'IBK_elderly', 
                                                             'LZ_elderly')] %>% 
    map(function(x) x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(gender$os$plots$IBK_young + 
                           theme(legend.position = 'right')), 
              ncol = 2, 
              rel_widths = c(0.9, 0.1))
  
  ## bottom panel
  
  paper_figures$gender_surv$bottom_panel <- list(gender$funct$plots$IBK_0$WHOFc + 
                                                   scale_y_continuous(breaks = 1:4, 
                                                                      labels = c('I', 'II', 'III', 'IV'), 
                                                                      limits = c(1, 4)), 
                                                 gender$funct$plots$LZ_0$WHOFc + 
                                                   scale_y_continuous(breaks = 1:4, 
                                                                      labels = c('I', 'II', 'III', 'IV'), 
                                                                      limits = c(1, 4)), 
                                                 gender$funct$plots$IBK_0$SMWD,
                                                 gender$funct$plots$LZ_0$SMWD, 
                                                 gender$funct$plots$IBK_0$NTproBNP + 
                                                   scale_y_continuous(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5), 
                                                                      trans = 'log10', 
                                                                      labels = c('10', '100', '1000', '10000', '100000'), 
                                                                      limits = c(1e1, 1e5)), 
                                                 gender$funct$plots$LZ_0$NTproBNP + 
                                                   scale_y_continuous(breaks = c(1e1, 1e2, 1e3, 1e4, 1e5), 
                                                                      trans = 'log10', 
                                                                      labels = c('10', '100', '1000', '10000', '100000'), 
                                                                      limits = c(1e1, 1e5))) %>% 
    map(function(x) x + theme(legend.position = 'none', 
                              axis.text.x = element_blank(), 
                              plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(gender$funct$plots$IBK_0$SMWD), 
              ncol = 2, 
              rel_widths = c(0.9, 0.1))
  
  ## complete figure
  
  paper_figures$gender_surv <- plot_grid(paper_figures$gender_surv$upper_panel, 
                                         paper_figures$gender_surv$bottom_panel, 
                                         nrow = 2, 
                                         rel_heights = c(0.4, 0.6), 
                                         labels = LETTERS, 
                                         label_size = 10) %>% 
    as_figure_object(figure_label = 'figure_5_gender_survival', 
                     w = 180, 
                     h = 230)
  
# Figure 6: Developed risk scores in the age and gender strata -----  
  
  insert_msg('Figure 6: risk scores in the age and gender strata')
  
  paper_figures$gender_signatures <- list(gender$funct$plots$IBK_0$sign_2525, 
                                          gender$funct$plots$LZ_0$sign_2525) %>% 
    map(function(x) x + theme(axis.text.x = element_blank(), 
                              legend.position = 'none', 
                              plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(gender$funct$plots$IBK_0$sign_2525), 
              ncol = 2, 
              rel_widths = c(0.9, 0.1)) %>% 
    as_figure_object(figure_label = 'figure_5_signatures_gender', 
                     w = 180, 
                     h = 80)
  
# Supplementary Figure S1: fit errors of the significant signatures in the train cohort, cv and test cohort -----
    
  insert_msg('Figure 5: redictiribution, cv and external validation errors')
  
  suppl_figures$fit_errors <- plot_grid(multi_plots$r_error_plots$development_mse + 
                                          labs(subtitle = 'Cox proportional hazard modeling'), 
                                        multi_plots$r_error_plots$testing_mse + 
                                          labs(subtitle = 'Cox proportional hazard modeling'), 
                                        ncol = 2, 
                                        align = 'hv', 
                                        labels = LETTERS, 
                                        label_size = 10) %>% 
    plot_grid(., multi_plots$c_plots$development_c + 
                labs(subtitle = 'Cox proportional hazard modeling, 20-fold CV'), 
              nrow = 2, 
              rel_heights = c(0.4, 0.6), 
              labels = c('', 'C'), 
              label_size = 10) %>% 
    as_figure_object(figure_label = 'figure_s1_model_errors', 
                     w = 180, 
                     h = 210)
  
# Supplementary Figure S2: Sensitivity, specificity and J statistic of the signatures at predicting 5-year mortality -----
  
  insert_msg('Figure S2: sensitivity, specificity and J statistic for the 5-year survival')
  
  suppl_figures$five_surv_stats <- plot_grid(multi_plots$five_surv_cutpoint$roc_sensitivity + 
                                               labs(subtitle = 'ROC') + 
                                               theme(legend.position = 'none'), 
                                             multi_plots$five_surv_cutpoint$roc_specificity + 
                                               labs(subtitle = 'ROC') + 
                                               theme(legend.position = 'none'), 
                                             multi_plots$five_surv_cutpoint$roc_j + 
                                               labs(subtitle = 'ROC') + 
                                               theme(legend.position = 'none'),
                                             get_legend(multi_plots$five_surv_cutpoint$roc_sensitivity), 
                                             ncol = 2, 
                                             align = 'hv', 
                                             axis = 'tblr', 
                                             labels = LETTERS, 
                                             label_size = 10) %>% 
    as_figure_object(figure_label = 'figure_s2_five_surv_prediction_stats', 
                     w = 180, 
                     h = 180)
  
# Supplementary Figure S3: Cox estimates, KM and ROC curves for the best-performing signature 2525 ----
  
  insert_msg('Figure S3: Cox estimates of the signatures of interest')
  
  ## top panel

  suppl_figures$best_sign$upper_panel <- plot_grid(multi_plots$signature_forests$sign_2525, 
                                                   ggdraw(), 
                                                   ncol = 2, 
                                                   align = 'hv', 
                                                   labels = LETTERS, 
                                                   label_size = 10)
  
  ## middle panel
  
  suppl_figures$best_sign$middle_panel <- list(km$plots$IBK_0$sign_2525, 
                                               km$plots$LZ_0$sign_2525) %>% 
  map(function(x) x + theme(legend.position = 'none', 
                            plot.tag = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(., get_legend(km$plots$IBK_0$sign_2525 + 
                              theme(legend.position = 'right')), 
              ncol = 2, 
              rel_widths = c(0.85, 0.15))
  
  ## bottom panel
  
  suppl_figures$best_sign$bottom_panel <- list(five_surv$roc_plots$IBK_0$sign_2525, 
                                               five_surv$roc_plots$LZ_0$sign_2525) %>% 
    map(function(x) x + theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(five_surv$roc_plots$IBK_0$sign_2525 + 
                           scale_color_manual(labels = c(globals$comp_labs, 
                                                         set_names('Sign 2525', 'sign_2525')), 
                                              values = c(globals$comp_colors, 
                                                         set_names('black', 'sign_2525')), 
                                              name = '')), 
              ncol = 2, 
              rel_widths = c(0.85, 0.15))
  
  ## complete figure
  
  suppl_figures$best_sign <- plot_grid(suppl_figures$best_sign$upper_panel, 
                                       suppl_figures$best_sign$middle_panel, 
                                       suppl_figures$best_sign$bottom_panel, 
                                       nrow = 3, 
                                       rel_heights = c(0.8, 1, 1), 
                                       labels = LETTERS, 
                                       label_size = 10) %>% 
    as_figure_object(figure_label = 'figure_s3_best_signature', 
                     w = 180, 
                     h = 210)

# Supplementary Figure S4: comparator scores in the age and gender strata -----
  
  insert_msg('Figure S4: comparator scores in the age and gender strata')
  
  suppl_figures$gender_comparators <- list(gender$funct$plots$IBK_0$mRASP, 
                                           gender$funct$plots$LZ_0$mRASP, 
                                           gender$funct$plots$IBK_0$Compera, 
                                           gender$funct$plots$LZ_0$Compera, 
                                           gender$funct$plots$IBK_0$SPAHR, 
                                           gender$funct$plots$LZ_0$SPAHR, 
                                           gender$funct$plots$IBK_0$FRENCH3p, 
                                           gender$funct$plots$LZ_0$FRENCH3p, 
                                           gender$funct$plots$IBK_0$FRENCH4p, 
                                           gender$funct$plots$LZ_0$FRENCH4p) %>% 
    map(function(x) x + theme(legend.position = 'none', 
                              plot.tag = element_blank(), 
                              axis.text.x = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(gender$funct$plots$IBK_0$mRASP), 
              ncol = 2, 
              rel_widths = c(0.9, 0.1)) %>% 
    as_figure_object(figure_label = 'figure_s4_comparator_gender', 
                     w = 180, 
                     h = 220)
  
# Supplementary Figure S5: Reveal scores in the age and gender strata -----
  
  insert_msg('Figure S5 Reveal scores in the age and gender strata')
  
  suppl_figures$gender_reveal <- list(gender$funct$plots$IBK_0$Reveal2_risk_3_cat, 
                                      gender$funct$plots$IBK_0$Reveal_lite2_3_cat) %>% 
    map(function(x) x + theme(legend.position = 'none', 
                              plot.tag = element_blank(), 
                              axis.text.x = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    plot_grid(get_legend(gender$funct$plots$IBK_0$Reveal2_risk_3_cat), 
              ncol = 2, 
              rel_widths = c(0.9, 0.1)) %>% 
    as_figure_object(figure_label = 'figure_s5_reveal_gender', 
                     w = 180, 
                     h = 60)
  
# saving the results on the disc ----
  
  insert_msg('Saving the plots')
  
  paper_figures %>% 
    walk(save_figure_object, 
         target_folder = './paper/figures', 
         device = cairo_pdf)

  suppl_figures %>% 
    walk(save_figure_object, 
         target_folder = './paper/supplementary figures', 
         device = cairo_pdf)
    
# END -----
  
  insert_tail()