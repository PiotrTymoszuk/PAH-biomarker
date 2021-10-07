# This script generates figures, tables and a report file

  c('./tools/sys_tools.R',
    './tools/plotting_tools.R', 
    './tools/project_tools.R') %>% 
    walk(source)
  
  library(writexl)
  
  insert_head()
  
# Rendering scripts -----
  
  insert_msg('Executing the paper scripts')
  
  c('./paper scripts/paper_figures.R',   
    './paper scripts/paper_tables.R', 
    './paper scripts/deploy_paper.R') %>% 
    walk(source)
  
# END ----
  
  insert_tail()