# This script reads and clears data from the PAH Innsbruck study

# libraries -----

  library(plyr)
  library(tidyverse)
  library(readxl)
  library(stringi)
  library(survival)
  library(soucer)
  library(figur)
  library(trafo)

  source_all('./tools/project_tools.R', 
             message = TRUE, 
             crash = TRUE)

  insert_head()

# data containers -----

  pah_study <- list()

# reading the data exported from SPSS, creating a legend data set. German typeset recoding ----
  
  insert_msg('Reading the raw data')

  pah_study$data_master <- 
    read_excel('./data/PAH_IbkandLinz_Wien_expforR_reveal.xlsx')
  
  ## there are issues with age at the first consultation
  
  pah_study$legend <- 
    read_excel('./data/Legend.xlsx')
  
# adding the center information to the patient ID and timepoint -----
  
  insert_msg('Updating ID and timepoint variables with center information')
  
  pah_study$data_master <- pah_study$data_master %>% 
    mutate(ID = ifelse(center == 1, paste('IBK', ID, sep = '_'), 
                       ifelse(center == 2, 
                              paste('LZ', ID, sep = '_'), 
                              paste('W', ID, sep = '_'))), 
           timepoint = ifelse(center == 1, 
                              paste('IBK', timepoint, sep = '_'), 
                              paste('LZ', timepoint, sep = '_')), 
           surv_months = Survival_time_from_FD_months)
  
# determining mortality within the study time frame, 1-, 3- and 5-year mortality -----
  
  insert_msg('Determining mortality within the study period')

  pah_study$data_master <- det_mortality(pah_study$data_master)
  
  pah_study$data_master <- pah_study$data_master %>% 
    mutate(event1 = factor(ifelse(death_study == 1 & surv_months < 12, 
                                  'yes', 'no')), 
           event3 = factor(ifelse(death_study == 1 & surv_months < 3 * 12, 
                                  'yes', 'no')), 
           event5 = factor(ifelse(death_study == 1 & surv_months < 5 * 12, 
                                  'yes', 'no')), 
           death_study_fct = factor(ifelse(death_study == 1, 
                                           'yes', 'no')))

# filtering out the CTEPH and the follow-ups, setting the variable format ----
  
  insert_msg('Formating and recoding the data')
  
  pah_study$data_master <- pah_study$data_master %>% 
    filter(PH_I_IV  == 0, 
           timepoint %in% c('IBK_0', 'LZ_0')) %>% 
    mutate(age_class = binarize(age_fc, 60, right = TRUE), 
           sex = car::recode(Gender, "0 = 'male'; 1 = 'female'") %>% 
             factor(c('female', 'male')), 
           mPAP_class = binarize(mPAP, 49, right = FALSE), ## https://doi.org/10.1093/ejechocard/jeq011
           dPAP_class = cut_median(dPAP, right = FALSE), ## https://doi.org/10.1093/ejechocard/jeq011
           TPG_class = cut_median(TPG, right = FALSE), 
           PVR_class = cut_median(PVR, right = TRUE),  ## https://journal.chestnet.org/article/S0012-3692(19)30152-7/pdf
           PCWP_class = binarize(PCWP, 12, right = TRUE),  ## https://www.ncbi.nlm.nih.gov/books/NBK557748/
           anemia = recode_01(anemia), 
           RDW_class = binarize(RDW, 14.5, right = TRUE), 
           RDW_log = log(RDW), 
           GFR_class = binarize(GFR, 60, right = FALSE) %>% 
             factor(c('\u226560', '<60')), 
           renal_ins = factor(ifelse(GFR < 60, 'yes', 'no'), 
                              c('no', 'yes')), 
           FT_class = cut(Ferritin, 
                          c(-Inf, 20, 250, Inf), 
                          c('\u226420', '20 - 250' , '>250')) %>% 
             factor(c('20 - 250', '\u226420', '>250')), 
           FT_log = log(Ferritin), 
           TSAT_class = binarize(TSAT, 20, right = FALSE) %>% 
             factor(c('<20', '\u226520')), 
           TSAT_log = log(TSAT), 
           iron_def_class = car::recode(iron_deficiency, 
                                        "0 = 'no'; 1 = 'mild'; 2 = 'severe'") %>% 
             factor(c('no', 'mild', 'severe')), 
           MCV_class = cut(MCV, 
                           c(-Inf, 80, 100, Inf), 
                           c('<80', '80 - 100' , '\u2265100'), 
                           right = F) %>% 
             factor(c('80 - 100', '<80', '\u2265100')), 
           NTproBNP_log = log(NTproBNP), 
           NTproBNP_class = cut(NTproBNP, ## https://journal.chestnet.org/article/S0012-3692(19)30152-7/pdf
                                c(-Inf, 300, 1100, Inf), 
                                c('<300', '300 - 1100' , '\u22651100')), 
           percardial_effusion = recode_01(percardial_effusion), 
           RAA_class = cut_median(RA_area, right = TRUE), 
           CI_class = cut(cardiac_index, ## ERS
                          c(-Inf, 2, 2.5, Inf), 
                          c('<2', '2 - 2.5', '\u22652.5'), 
                          right = F), 
           mRAP_class = cut(mRAP, ## ERS
                            c(-Inf, 8, 14, Inf), 
                            c('<8', '8 - 14', '>14')), ## ERS
           WHOFc_class = cut(WHOFc, 
                             c(-Inf, 2.5, Inf), 
                             c('I/II', 'III/IV')), ## https://journal.chestnet.org/article/S0012-3692(19)30152-7/pdf
           SMWD_class = cut(SMWD, 
                            c(-Inf, 165, 320, 440, Inf), 
                            c('<165', '165 - 320', '321 - 439', '\u2265440'), 
                            right = F) %>% 
             factor(c('\u2265440', '321 - 439', '165 - 320', '<165')), 
           SO2_RL_class = binarize(SO2_RL, 95, right = FALSE) %>% 
             factor(c('\u226595', '<95')), 
           SvO2_class = cut(SvO2, 
                            c(-Inf, 60, 65, Inf), 
                            c('<60', '60 - 65', '>65'), 
                            right = F) %>% 
             factor(c('>65', '60 - 65', '<60')), ## ERS
           CRP_class = binarize(CRP, 0.5, right = F), 
           HCT_class = binarize(HCT, 0.40, right = F) %>% 
             factor(c('<0.4', '\u22650.4')), 
           creatinine_class = ifelse(sex == 'male', 
                                     ifelse(Creatinin > 1.3, 'high', 'normal'), 
                                     ifelse(Creatinin > 1.1, 'high', 'normal')) %>% 
             factor(c('normal', 'high')), 
           pO2_RL_class = binarize(pO2_RL, 80, right = FALSE) %>% 
             factor(c('\u226580', '<80')), 
           DLCO_class = binarize(DLCO, 80, right = FALSE) %>% 
             factor(c('\u226580', '<80')), 
           KCO_class = binarize(KCO, 80, right = FALSE) %>% 
             factor(c('\u226580', '<80')), 
           TAPSE_class = binarize(TAPSE, 17, right = FALSE) %>% 
             factor(c('\u226517', '<17')), 
           LVEF_class = binarize(LVEF, 70, right = FALSE) %>% 
             factor(c('\u226570', '<70')), 
           BMI_class = cut(BMI, 
                           c(-Inf, 25, 30, Inf), 
                           c('normal', 'overweight', 'obesity')))
  
# generating the vectors with modeling variables (present in both cohorts) ----
    
  insert_msg('Selecting the modeling responses and variables')

  ## responses, independent variables and comparators
  
  pah_study[c('mod_variables', 
              'responses', 
              'comparators')] <- c('independent', 
                                   'response', 
                                   'comparator') %>% 
    map(~filter(pah_study$legend, variable_type == .x))
  
  ## non-baseline variable levels 
  
  pah_study$levels <- pah_study$mod_variables$variable %>% 
    map(~levels(pah_study$data_master[[.x]])) %>% 
    map(~.x[2:length(.x)]) %>% 
    set_names(pah_study$mod_variables$variable)

# recoding risk scales -------
  
  insert_msg('Recoding risk scales')
  
  pah_study$data_master <- pah_study$data_master %>% 
    mutate(mRASP = car::recode(mRASP, 
                               "0 = 'low'; 1 = 'int'; 2 = 'high'"), 
           mRASP = factor(mRASP, 
                          c('low', 'int', 'high')), 
           Compera = car::recode(Compera, 
                                 "1 = 'low'; 2 = 'int'; 3 = 'high'"), 
           Compera = factor(Compera, 
                            c('low', 'int', 'high')), 
           SPAHR  = car::recode(SPAHR , 
                                "1 = 'low'; 2 = 'int'; 3 = 'high'"), 
           SPAHR  = factor(SPAHR , 
                           c('low', 'int', 'high')), 
           FRENCH3p = factor(FRENCH3p), 
           FRENCH4p = factor(FRENCH4p), 
           Reveal_lite2_3_cat  = car::recode(Reveal_lite2_3_cat , 
                                             "1 = 'low'; 2 = 'int'; 3 = 'high'"), 
           Reveal_lite2_3_cat  = factor(Reveal_lite2_3_cat , 
                                        c('low', 'int', 'high')), 
           Reveal2_risk_3_cat  = car::recode(Reveal2_risk_3_cat , 
                                             "1 = 'low'; 2 = 'int'; 3 = 'high'"), 
           Reveal2_risk_3_cat  = factor(Reveal2_risk_3_cat , 
                                        c('low', 'int', 'high')))
  
# dividing the dataset into a list of datasets for to the particular time points ----

  insert_msg('Creating data sheets for each center')
  
  pah_study <- pah_study$data_master %>% 
    blast(timepoint) %>% 
    c(pah_study, .)

# some globals: color scales and plot labels ----
  
  insert_msg('Creating a list with plotting and modeling globals')
  
  source_all('./tools/project_globals.R', 
             message = TRUE, crash = TRUE)
  
# generating variable: level dictionary -----
  
  pah_study$level_dict <- pah_study$data_master %>% 
    map(function(x) if(is.factor(x)) levels(x) else NULL) %>% 
    compact %>% 
    map2_dfr(., names(.), ~tibble(variable = .y, level = .x)) %>% 
    mutate(parameter = paste0(variable, level),
           parameter = stri_replace(parameter, 
                                    regex = '\u2264|\u2265', 
                                    replacement = '='), 
           var_label = exchange(variable, 
                                dict = globals$var_labs), 
           param_label = paste(var_label, level, sep = ': '))
  
# Data sets for IBK and LZ, first consultation with the modeling variables and complete cases only ----
  
  insert_msg('Modeling data sets')
  
  pah_study[c('IBK_0', 'LZ_0')] <- pah_study[c('IBK_0', 'LZ_0')] %>% 
    map(select, 
        ID, 
        event3, 
        event5, 
        death_study, 
        death_study_fct, 
        surv_months, 
        all_of(pah_study$mod_variables$variable)) %>% 
    map(~filter(.x, complete.cases(.x)))
  
# creating a series of relevant survival objects: IBK timepoints 0 and 1 and Linz timepoint 0 ----
  
  insert_msg('Survival objects for the relevant timepoints and centers')
  
  pah_study$surv_obj <- pah_study[c('IBK_0', 
                                    'LZ_0')] %>% 
    map(~Surv(time = .x[['surv_months']], 
              event = .x[['death_study']]))
  
# END ----
  
  insert_tail()