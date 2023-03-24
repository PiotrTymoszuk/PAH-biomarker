# Performs data analysis

# tools ----

  library(soucer)
  library(rms)
  library(rlang)
  library(glmnet)
  library(caret)
  library(furrr)
  library(survminer)
  library(coxExtensions)
  library(exda)
  library(rstatix)
  library(DescTools)
  library(somKernels)
  library(clustTools)
  library(ggrepel)
  library(psych)
  library(randomForestSRC)
  library(pec)

  c('./tools/project_tools.R', 
    './tools/project_globals.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
  insert_head()
  
# Analysis scripts -----
  
  insert_msg('Analysis scripts')
  
  source_all(c('./analysis scripts/pca.R', 
               ## univariate models and Elastic Net modeling
               './analysis scripts/univariate_cox.R', 
               './analysis scripts/multivariate_cox.R', 
               './analysis scripts/multi_visualization.R', 
               ## performance and ensemble of risk assessment tools
               './analysis scripts/tool_correlation.R',
               './analysis scripts/tool_ensemble.R', 
               './analysis scripts/tool_lasso.R', 
               './analysis scripts/tool_comparison.R', 
               ## clustering
               './analysis scripts/cluster_devel.R', 
               './analysis scripts/clustering.R', 
               './analysis scripts/cluster_characteristic.R', 
               './analysis scripts/cluster_survival.R', 
               './analysis scripts/cluster_reviewer.R'), 
             message = TRUE, crash = TRUE)

# END ----
  
  insert_tail()