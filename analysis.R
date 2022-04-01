# Performs data analysis -----

  library(soucer)
  library(rms)
  library(rlang)
  library(glmnet)
  library(caret)
  library(furrr)
  library(survminer)
  library(coxExtensions)
  library(exda)
  library(somKernels)
  library(clustTools)
  library(ggrepel)
  
  source_all(c('./tools/project_tools.R'), 
             message = TRUE, crash = TRUE)
  
  insert_head()
  
# Analysis scripts -----
  
  source_all(c('./analysis scripts/pca.R', 
               './analysis scripts/univariate_cox.R', 
               './analysis scripts/multivariate_cox.R', 
               './analysis scripts/multi_visualization.R', 
               './analysis scripts/cluster_devel.R', 
               './analysis scripts/clustering.R', 
               './analysis scripts/cluster_characteristic.R', 
               './analysis scripts/cluster_survival.R'), 
             message = TRUE, crash = TRUE)

# END ----
  
  insert_tail()