# this script renders a word file with the figures and supplementary material

  insert_head()
  
# Exporting the figure Rmd chunks ----
  
  insert_msg('Figure chunks')

  insert_figure(paper_figures$strobe, 
                paper_figures$elanet, 
                paper_figures$part_clust, 
                paper_figures$cluster_risk, 
                paper_figures$summary, 
                file = './paper/markdown/figure_chunks.Rmd', 
                ref_names = stri_replace(names(paper_figures), fixed = '_', replacement = '-'), 
                captions = c('CONSORT flow diagram of the study analysis inclusion process.', 
                             'Multi-parameter survival modeling.', 
                             'Clustering of the study participants.', 
                             'Risk assessment and survival differences in the participant clusters.', 
                             'Summary of the analysis results.'), 
                add_extern_legend = TRUE, 
                append = FALSE)
    
# exporting the supplementary figure Rmd chunks ------
  
  insert_msg('Supplementary Figure chunks')
  
  insert_figure(suppl_figures$uni_cox, 
                suppl_figures$elanet_calibration, 
                suppl_figures$cluster_qc, 
                suppl_figures$clust_diff, 
                file = './paper/markdown/suppl_chunks.Rmd', 
                ref_names = stri_replace_all(names(suppl_figures), fixed = '_', replacement = '-'), 
                captions = c('Univariable Cox proportional hazard modeling.', 
                             'Elastic net model linear prediction score and overall survival.', 
                             'Development of participant clusters.', 
                             'Differences in study variable between the participant clusters.'), 
                add_extern_legend = TRUE, 
                append = FALSE)
  
# rendering the figures and tables ------

  insert_msg('Rendering the figures and tables')
  
  render('./paper/markdown/figures_and_tables.Rmd', 
         output_format = word_document2(number_sections = FALSE, 
                                        reference_docx = 'ms_template.docx'), 
         output_dir = './paper/') 

# rendering the supplementary material -----
  
  insert_msg('Rendering the supplementary material')
  
  render('./paper/markdown/supplementary_material.Rmd', 
         output_format = word_document2(number_sections = FALSE, 
                                        reference_docx = 'ms_template.docx'), 
         output_dir = './paper/') 
  
# END ------
  
  insert_tail()