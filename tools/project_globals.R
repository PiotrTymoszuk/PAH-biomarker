# The project globals

  insert_head()
  
# container ------
  
  globals <- list()
  
# plotting theme ------
  
  globals$common_text <- 
    element_text(size = 8, face = 'plain', color = 'black')
  
  globals$common_margin <- 
    ggplot2::margin(t = 4, l = 3, r = 2, unit = 'mm')
  
  globals$common_theme <- theme_classic() + 
    theme(axis.text = globals$common_text, 
          axis.title = globals$common_text, 
          plot.title = element_text(size = 8, 
                                    face = 'bold', 
                                    color = 'black', 
                                    hjust = 0), 
          plot.subtitle = globals$common_text, 
          plot.tag = element_text(size = 8, 
                                  face = 'plain', 
                                  color = 'black', 
                                  hjust = 0), 
          plot.tag.position = 'bottom', 
          legend.text = globals$common_text, 
          legend.title = globals$common_text, 
          strip.text = globals$common_text,
          strip.background = element_rect(fill = 'gray95', color = 'gray80'), 
          plot.margin = globals$common_margin, 
          panel.grid.major = element_line(color = 'gray90'))
  
# plotting colors ------
  
  globals$cluster_colors <- c(`#1` = 'coral3', 
                              `#2` = 'steelblue', 
                              `#3` = 'gray60', 
                              `#4` = 'darkolivegreen3')
  
  globals$pos_neg_scale <- c('negative' = 'steelblue', 
                             'positive' = 'coral3', 
                             'ns' = 'gray60')
  
  globals$center_colors <- c(IBK_0 = 'coral2', 
                             cv = 'darkorange4', 
                             LZ_0 = 'lightskyblue3')
  
  globals$center_labs <- c(IBK_0 = 'IBK', 
                           cv = 'CV IBK', 
                           LZ_0 = 'LZ/W')
  
# variable lexicon --------
  
  globals$var_labs <- pah_study$mod_variables$label %>% 
    set_names(pah_study$mod_variables$variable) %>% 
    c(., c(Gendermale = 'Sex: male', 
           percardial_effusionyes = 'Percardial\neffusion', 
           mRASP = 'mRASP', 
           Compera = 'Compera', 
           SPAHR = 'SPAHR', 
           FRENCH3p = 'FRENCH3p', 
           FRENCH4p = 'FRENCH4p', 
           Reveal_lite2_3_cat = 'Reveal Lite', 
           Reveal2_risk_3_cat = 'Reveal 2.0')) %>% 
    compress(names_to = 'variable', 
             values_to = 'label')
  
# END ------
  
  insert_tail()
