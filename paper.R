# Generates paper tables and figures

  library(knitr)
  library(rmarkdown)
  library(bookdown)
  library(flextable)
  library(writexl)
  
  library(soucer)
  library(figur)
  
  source_all(c('./tools/project_tools.R'), 
             message = TRUE, crash = TRUE)
  
  insert_head()
  
# paper scripts -----
  
  insert_msg('Sourcing the paper scripts')
  
  source_all(c('./paper scripts/paper_tables.R', 
               './paper scripts/paper_figures.R', 
               './paper scripts/deploy_paper.R'), 
             message = TRUE, crash = TRUE)
  
# END -----
  
  insert_tail()