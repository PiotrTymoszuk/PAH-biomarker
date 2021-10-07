# This script reads and clears data from the PAH Innsbruck study

# libraries -----

  library(plyr)
  library(tidyverse)
  library(readxl)
  library(stringi)
  library(survival)

  c('./tools/sys_tools.R', 
    './tools/km_toolbox.R', 
    './tools/plotting_tools.R', 
    './tools/project_tools.R') %>% 
    walk(source)

  insert_head()

# data containers -----

  pah_study <- list()
  globals <- list()

# reading the data exported from SPSS, creating a legend data set. German typeset recoding ----
  
  insert_msg('Reading the raw data')

  pah_study$data_master <- read_excel('./input data/PAH_IbkandLinz_Wien_expforR_reveal.xlsx')
  
  pah_study$legend <- read_excel('./input data/Legend.xlsx')
  
# adding the center information to the patient ID and timepoint -----
  
  insert_msg('Updating patient ID and timepoint variables with center information')
  
  pah_study$data_master <- pah_study$data_master %>% 
    mutate(ID = ifelse(center == 1, paste('IBK', ID, sep = '_'), 
                       ifelse(center == 2, paste('LZ', ID, sep = '_'), paste('W', ID, sep = '_'))), 
           timepoint = ifelse(center == 1, paste('IBK', timepoint, sep = '_'), paste('LZ', timepoint, sep = '_')))
  
# determining mortality within the study time frame and 5-year mortality -----
  
  insert_msg('Determining mortality within the study period')

  pah_study$data_master <- det_mortality(pah_study$data_master)
  
  pah_study$data_master <- pah_study$data_master %>% 
    mutate(death_acute = ifelse(death_study == 1 & Survival_time_from_FD_months < 5 * 12, 1, 0))
  
# filtering out the CTEPH and the follow-ups, setting the variable format ----
  
  insert_msg('Formating and recoding the data')
  
  pah_study$data_master <- pah_study$data_master %>% 
    filter(PH_I_IV  == 0, 
           timepoint %in% c('IBK_0', 'LZ_0')) %>% 
    mutate(age_class = binarize(age_fc, 65, right = T), 
           sex = car::recode(Gender, "0 = 'male'; 1 = 'female'") %>% 
             factor(c('male', 'female')), 
           mPAP_class = binarize(mPAP, 49, right = F), ## https://doi.org/10.1093/ejechocard/jeq011
           dPAP_class = cut_median(dPAP, right = F), ## https://doi.org/10.1093/ejechocard/jeq011
           TPG_class = cut_median(TPG, right = F), 
           PVR_class = cut_median(PVR, right = T),  ## https://journal.chestnet.org/article/S0012-3692(19)30152-7/pdf
           PCWP_class = binarize(PCWP, 12, right = T),  ## https://www.ncbi.nlm.nih.gov/books/NBK557748/
           anemia = recode_01(anemia), 
           RDW_class = binarize(RDW, 14.5, right = T), 
           GFR_class = binarize(GFR, 60, right = F) %>% 
             factor(c('\u226560', '<60')), 
           FT_class = cut(Ferritin, 
                          c(-Inf, 20, 250, Inf), 
                          c('\u226420', '20 - 250' , '>250')) %>% 
             factor(c('20 - 250', '\u226420', '>250')), 
           TSAT_class = binarize(TSAT, 20, right = F) %>% 
             factor(c('<20', '\u226520')), 
           iron_def_class = car::recode(iron_deficiency, 
                                        "0 = 'no'; 1 = 'mild'; 2 = 'severe'") %>% 
             factor(c('no', 'mild', 'severe')), 
           MCV_class = cut(MCV, 
                           c(-Inf, 80, 100, Inf), 
                           c('<80', '80 - 100' , '\u2265100'), 
                           right = F) %>% 
             factor(c('80 - 100', '<80', '\u2265100')), 
           NTproBNP_class = cut(NTproBNP, ## https://journal.chestnet.org/article/S0012-3692(19)30152-7/pdf
                                c(-Inf, 300, 1100, Inf), 
                                c('<300', '300 - 1100' , '\u22651100')), 
           percardial_effusion = recode_01(percardial_effusion), 
           RAA_class = cut_median(RA_area, right = T), 
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
           SO2_RL_class = binarize(SO2_RL, 95, right = F) %>% 
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
           pO2_RL_class = binarize(pO2_RL, 80, right = F) %>% 
             factor(c('\u226580', '<80')), 
           DLCO_class = binarize(DLCO, 80, right = F) %>% 
             factor(c('\u226580', '<80')), 
           KCO_class = binarize(KCO, 80, right = F) %>% 
             factor(c('\u226580', '<80')), 
           TAPSE_class = binarize(TAPSE, 17, right = F) %>% 
             factor(c('\u226517', '<17')), 
           LVEF_class = binarize(LVEF, 70, right = F) %>% 
             factor(c('\u226570', '<70')), 
           BMI_class = cut(BMI, 
                           c(-Inf, 25, 30, Inf), 
                           c('normal', 'overweight', 'obesity')))
  
# generating the vectors with modeling variables (present in both cohorts) ----
    
  insert_msg('Selecting the modeling responses and variables')
  
  ## variables of interest
  
  pah_study$data_master <- pah_study$data_master %>% 
    select(all_of(pah_study$legend$variable))
  
  ## responses, independent variables and comparators
  
  pah_study[c('mod_variables', 
              'responses', 
              'comparators')] <- c('independent', 
                                   'response', 
                                   'comparator') %>% 
    map(function(x) filter(pah_study$legend, variable_type == x))
  
  ## non-baseline variable levels 
  
  pah_study$levels <- pah_study$mod_variables$variable %>% 
    map(function(x) levels(pah_study$data_master[[x]])) %>% 
    map(function(x) x[2:length(x)]) %>% 
    set_names(pah_study$mod_variables$variable)

# dividing the dataset into a list of datasets corresponding to the particular time points ----

  insert_msg('Creating data sheets for each center')
  
  pah_study <- pah_study$data_master %>% 
    dlply(.(timepoint), as_tibble) %>% 
    c(pah_study, .)
  
# creating a series of relevant survival objects: IBK timepoints 0 and 1 and Linz timepoint 0 ----
  
  insert_msg('Generating a series of survival objects for the relevant timepoints and centers')
  
  pah_study$surv_obj <- pah_study[c('IBK_0', 
                                    'LZ_0')] %>% 
    map(create_surv, 
        time_variable = 'Survival_time_from_FD_months', 
        event_variable = 'death_study')
  
# a table with centers and their plotting order -----
  
  insert_msg('Preparing a table with center information for plotting')
  
  pah_study$center_tbl <- tibble(center = c('IBK_0', 'LZ_0'), 
                                 center_lab = c('IBK', 'LZ/W'), 
                                 center_order = c(1, 2))

# some globals: color scales and plot labels ----
  
  insert_msg('Creating a list with plotting and modeling globals')
  
  globals$pos_neg_scale <- c('negative' = 'steelblue', 
                             'positive' = 'coral3', 
                             'ns' = 'gray60')

  globals$center_colors <- c(IBK_0 = 'coral2', 
                             cv = 'darkorange4', 
                             LZ_0 = 'lightskyblue3')
  
  globals$center_labs <- c(IBK_0 = 'IBK', 
                           cv = 'CV IBK', 
                           LZ_0 = 'LZ/W')
  
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
           Reveal2_risk_3_cat = 'Reveal 2.0'))

  globals$mort_colors <- c('0' = 'cornsilk', 
                           '1' = 'firebrick4')
  
  globals$mort_labs <- c('0' = 'Survived', 
                         '1' = 'Deceased')
  
  
  ## theme
  
  globals$common_text <- element_text(size = 8, face = 'plain', color = 'black')
  
  globals$common_margin <- ggplot2::margin(t = 4, l = 3, r = 2, unit = 'mm')
  
  globals$common_theme <- theme_classic() + theme(axis.text = globals$common_text, 
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
                                                  plot.margin = globals$common_margin)
  
  ## comparator labels and colors
  
  globals$comp_colors <- c('darkseagreen4', 
                           'steelblue3', 
                           'steelblue4', 
                           'plum2', 
                           'plum4', 
                           'darkgoldenrod2', 
                           'darkgoldenrod4') %>% 
    set_names(pah_study$comparators$variable)
  
  globals$comp_labs <- set_names(pah_study$comparators$label, 
                                 pah_study$comparators$variable)
  
  ## signature type color and labels
  
  globals$signature_color <- c(new = 'steelblue', 
                               comparator = 'coral3')
  
  globals$signature_labels <- c(new = 'Developed', 
                                comparator = 'Comparator')
  
# generating variable: level dictionary -----
  
  pah_study$level_dict <- pah_study$data_master %>% 
    map(function(x) if(is.factor(x)) levels(x) else NULL) %>% 
    compact %>% 
    map2_dfr(., names(.), ~tibble(variable = .y, level = .x)) %>% 
    mutate(parameter = paste0(variable, level),
           parameter = stri_replace(parameter, regex = '\u2264|\u2265', replacement = '='), 
           var_label = translate_vars(variable), 
           param_label = paste(var_label, level, sep = ': '))
  
# END ----
  
  insert_tail()