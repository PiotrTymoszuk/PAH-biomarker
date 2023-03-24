# Links to development packages ------

  insert_head()

# container -----

  links <- list()
  
# links -----
  
  insert_msg('Mdlink objects')
  
  links <- c('_trafo_' = 'https://github.com/PiotrTymoszuk/trafo', 
             '_figur_' = 'https://github.com/PiotrTymoszuk/figur', 
             '_coxExtensions_' = 'https://github.com/PiotrTymoszuk/coxExtensions', 
             '_ExDA_' = 'https://github.com/PiotrTymoszuk/ExDA', 
             '_clustTools_' = 'https://github.com/PiotrTymoszuk/clustTools') %>% 
    compress(names_to = 'ref_name', 
             values_to = 'x') %>%
    pmap(mdlink) %>% 
    set_names(c('trafo', 'figur', 'coxExtensions', 'exda', 'clustTools'))
    
  
# END ------
  
  insert_tail()