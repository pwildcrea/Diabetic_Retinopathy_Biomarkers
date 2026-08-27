##### - NDD Blood-Based Biomarkers & Diabetic Retinopathy Assessed by ETDRS Data Analysis - #####
install.packages("dplyr")
install.packages("corrr")
install.packages("ggcorrplot")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("ggplot2")
install.packages("haven")
install.packages("ggpubr")
install.packages("devtools")
install.packages("ggeffects")
install.packages("modelsummary")
install.packages("gtsummary")
install.packages("lmtest")
install.packages("quantreg")
install.packages("table1")
install.packages("Hmisc")
install.packages("haven")
install.packages("tidyverse")
install.packages("forestplot")
install.packages("bestNormalize")
install.packages("ordinal")
library(dplyr)
library(corrr)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(haven)
library(ggpubr)
library(ggeffects)
library(modelsummary)
library(gtsummary)
library(lmtest)
library(quantreg)
library(table1)
library(Hmisc)
library(haven)
library(tidyverse)
library(forestplot)
library(broom)
library(grid)
library(bestNormalize)
library(ordinal)

#######################################################################################################

# Data Cleaning & Structuring

#######################################################################################################
# Load Files
baseline_characteristics <- read.csv("/Users/Peter/Desktop/UMN Medical School/Research/R21/baseline_characteristics.csv", header = TRUE)
baseline_MRI <- read.csv("/Users/Peter/Desktop/UMN Medical School/Research/R21/baseline_MRI.csv", header = TRUE)
all_MRI <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/Mind_MRI_FINAL.csv')
ev_id <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/linking_id_ev.csv')
ev_prot <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/ndx_data_clean.csv')
a1c_accordion <- read_sas('/Users/Peter/Desktop/UMN Medical School/Research/R21/ACCORD_2017b_2 (2)/Ancillary_Study/ACCORDION/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/hba1c.sas7bdat')
baseline_hx <- read_sas('/Users/Peter/Desktop/UMN Medical School/Research/R21/ACCORD_2017b_2 (2)/Main_Study/4-Data_Sets-CRFs/4a-CRF_Data_Sets/f07_baselinehistoryphysicalexam.sas7bdat')
a1c_subset <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/a1c_subset.csv')
accord_race <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/ACCORD_2017b_2 (2)/Main_Study/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/csv/accord_key.csv')
#f80_a1c <-  read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/F80_hba1c.csv')
detailed_mri <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/MUSE_MRI_DERIVED_VOLUMES.csv')
detailed_mri_ids <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/detailed_mri_ids.csv')
mri_detailed_blr <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/detailed_mri_blr.csv')
mri_detailed_f40 <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/detailed_mri_f40.csv')
mri_detailed_f80 <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/detailed_mri_f80.csv')
detailed_mri_ids <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/detailed_mri_ids.csv')
stroke_data <- read_sas('/Users/Peter/Desktop/UMN Medical School/Research/R21/ACCORD_2017b_2 (2)/Main_Study/4-Data_Sets-CRFs/4a-CRF_Data_Sets/f01_inclusionexclusionsummary.sas7bdat')
# Merge Only Baseline Data
baseline_C_MRI <- merge(baseline_characteristics, baseline_MRI,
                        by.x = "MaskID")
baseline_merge <- merge(baseline_C_MRI, baseline_hx,
                        by.x = "MaskID")
# Calculate BMI as New Column
baseline_merge$BMI <- (baseline_merge$wt_kg / (baseline_merge$ht_cm/100)^2)
# Merged Characteristics, All MRI Data, EV IDs, & EV Protein
#ev_id = rename(ev_id, "prot_ID" = "prot_id")                         #Rename ID column 
ev_prot <- merge(ev_id, ev_prot,
                 by.x = "prot_ID")                             #Merge by IDs
ev_prot <- ev_prot[, -which(names(ev_prot) == "prot_ID")]
ev_prot <- ev_prot[, -which(names(ev_prot) == "sample")]             #Delete additional ID column
# Replace NaN w/ 0
ev_prot[is.na(ev_prot)] <- 0
# Generate Parent File With All Data
baseline_all <- merge(ev_prot, baseline_merge,
                      by.x = "MaskID")
# Convert pTau Data to Lumipulse Equivalent Values
baseline_all$lmp_ptau_217 <- ((0.0026 * baseline_all$ptau217 ) + 0.0586)
# Subset Accordion A1c Data
tau_to_merge <- tibble(baseline_all$MaskID, baseline_all$log_lmp_ptau_217)
a1c_accordion = rename(a1c_accordion, "MaskID" = "MASKID")  
a1c_accordion_subset <- merge(a1c_accordion, ev_id,
                              by.x = "MaskID")
tau_to_merge = rename(tau_to_merge, "MaskID" = "baseline_all$MaskID")  
a1c_subset <- merge(a1c_subset, tau_to_merge,
                    by.x = "MaskID")
# Add ACCORD race data to baseline_all
accord_race <- accord_race[,c("MaskID", "raceclass")]
baseline_all <- merge(baseline_all, accord_race,
                      by.x = "MaskID")
# Switch ACCORDION "race" to "raceclass" and ACCORD "raceclass" to "race"
baseline_all = rename(baseline_all, "racee" = "raceclass")
baseline_all = rename(baseline_all, "raceclass" = "race")
baseline_all = rename(baseline_all, "race" = "racee")
# Convert Tx to Binary Standard vs Intensive
baseline_all <- baseline_all %>%
  mutate(treatment = (recode(treatment, 
                             "Intensive Gylcemia/Standard BP" = "Intensive Glycemic Control",
                             "Intensive Gylcemia/Intensive BP" = "Intensive Glycemic Control",
                             "Intensive Glycemia/Lipid Placebo" = "Intensive Glycemic Control",
                             "Intensive Glycemia/Lipid Fibrate" = "Intensive Glycemic Control",
                             "Standard Gylcemia/Standard BP" = "Standard Glycemic Control",
                             "Standard Gylcemia/Intensive BP" = "Standard Glycemic Control",
                             "Standard Glycemia/Lipid Placebo" = "Standard Glycemic Control",
                             "Standard Glycemia/Lipid Fibrate" = "Standard Glycemic Control",)))
# Replace pTau values below LLOQ w/ Protocol LLOQ 1.95
ptau_LLQ_threshold <- 1.95
baseline_all$ptau217 <- replace(baseline_all$ptau217, baseline_all$ptau217< 1.95, 1.95)
baseline_all$ptau217[is.na(baseline_all$ptau217)] <- 1.95
# Transform Protein & Yrs Diabetes Data w/ Natural Log (C + 1)
baseline_all$logAb40 <- log(baseline_all$Ab40 + 1)
baseline_all$logAb42 <- log(baseline_all$Ab42 + 1)
baseline_all$logAb38 <- log(baseline_all$Ab38 + 1)
baseline_all$logGSK3B <- log(baseline_all$GSK3B + 1)
baseline_all$logIGF1R <- log(baseline_all$IGF1R + 1)
baseline_all$logIRS1 <- log(baseline_all$IRS1 + 1)
baseline_all$logAKT <- log(baseline_all$AKT + 1)
baseline_all$logmTOR <- log(baseline_all$mTOR + 1)
baseline_all$logp70S6K <- log(baseline_all$p70S6K + 1)
baseline_all$logIR <- log(baseline_all$IR + 1)
baseline_all$logPTEN <- log(baseline_all$PTEN + 1)
baseline_all$logGSK3a <- log(baseline_all$GSK3a + 1)
baseline_all$logTSC2 <- log(baseline_all$TSC2 + 1)
baseline_all$logRPS6 <- log(baseline_all$RPS6 + 1)
baseline_all$logptau217 <- log(baseline_all$ptau217 + 1)
baseline_all$log_lmp_ptau_217 <- log(baseline_all$lmp_ptau_217 + 1)
baseline_all$logNFL <- log(baseline_all$NFL + 1)
baseline_all$logyrsdiab <- log(baseline_all$yrsdiab + 1)
# Compute Ab4240 Ratio w/ Log Transformed Data
baseline_all$logAb4240_Ratio <- (baseline_all$logAb42 / baseline_all$logAb40)
# Compute pTau-217/Ab4240 Ratio w/ Log Transformed Data
# Arbitrary Units
baseline_all$logpT217_AB42_Ratio <- (baseline_all$logptau217 / baseline_all$logAb42)
# pg/dL
baseline_all$logLMP_pT217_AB42_Ratio <- (baseline_all$log_lmp_ptau_217 / baseline_all$logAb42)
# Subset Follow Up Data
mind_base <- all_MRI[all_MRI$VISIT == 'BLR',]
mind_F40 <- all_MRI[all_MRI$VISIT == 'F40',]
mind_F80 <- all_MRI[all_MRI$VISIT == 'F80',]
# Add Visit Suffix To All Column Names
colnames(mind_base) <- paste(colnames(mind_base), "BLR", sep = "_")
colnames(mind_F40) <- paste(colnames(mind_F40), "F40", sep = "_")
colnames(mind_F80) <- paste(colnames(mind_F80), "F80", sep = "_")
View(mind_F80)        #Just to check addition of suffix
# Rename ID column
mind_base = rename(mind_base, "MaskID" = "ID_BLR")
mind_F40 = rename(mind_F40, "MaskID" = "ID_F40")
mind_F80 = rename(mind_F80, "MaskID" = "ID_F80")
View(mind_F80)       # Just to check change of ID name
# Remove Match column
mind_base <- subset(mind_base, select = -ID_Match_BLR)
mind_F40 <- subset(mind_F40, select = -ID_Match_F40)
mind_F80 <- subset(mind_F80, select = -ID_Match_F80)
# Reload correct files after editing in excel and removing duplicates/transfering from F40 to F80
baseline_all_blr <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/baseline_all_blr.csv')
baseline_all_f40 <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/baseline_all_f40.csv')
baseline_all_f80 <- read.csv('/Users/Peter/Desktop/UMN Medical School/Research/R21/baseline_all_f80.csv')
# Merge baseline_all and all 3 MRI follow up datasets
baseline_all <- baseline_all %>%
  left_join(baseline_all_blr, by = "MaskID") %>%
  left_join(baseline_all_f40, by = "MaskID") %>%
  left_join(baseline_all_f80, by = "MaskID") %>%
  left_join(mind_base, by = "MaskID") %>%
  left_join(mind_F40, by = "MaskID") %>%
  left_join(mind_F80, by = "MaskID")
View(baseline_all)
# Relevel Race for MLR
baseline_all |> count(race)
# Make factor
baseline_all$race <- factor(baseline_all$race)
# Relevel
baseline_all$race <- relevel(baseline_all$race, ref = "White")
View(baseline_all)
# Subset Stroke Data
stroke <- stroke_data[c("MaskID", "x2stroke")]
View(stroke)
# Add Stroke Data to baseline_all
baseline_all <- merge(baseline_all, stroke,
                    by.x = "MaskID")
# Make Any Diabetes Medication Use a Column
baseline_all <- baseline_all %>%
  mutate(
    any_diab_med = as.integer(
      rowSums(across(c(88:99)), na.rm = TRUE) > 0
    )
  )

#######################################################################################################

# Subset Data to Produce Cleaner Data Frame

#######################################################################################################
# Subset Data
retinsadf <- baseline_all %>% 
  select("MaskID", 
         "baseline_age", 
         "sex", 
         "edu", 
         "race", 
         "cigarett",
         "yrsdiab",
         "linking_hba1cAccordion.HBA1C", 
         "Linking_blodpressure.Accordion.SBP", 
         "BMI", 
         "Linking_Lipids.Accordion.CHOL", 
         "cvd_hx_baseline", 
         "Linking.otherlabs.Accordion.GFR",
         "Linking_Lipids.Accordion.LDL",
         "Linking_Lipids.Accordion.HDL",
         "Linking_blodpressure.Accordion.DBP",
         "alcohol",
         "Linking_concomitantmedsAccordio.MEM_LOSS",
         "Linking_concomitantmedsAccordio.ANTI_DEPRESS",
         "ptau217", 
         "lmp_ptau_217", 
         "NFL",
         "logptau217", 
         "log_lmp_ptau_217", 
         "logNFL", 
         "logAb4240_Ratio", 
         "logpT217_AB42_Ratio", 
         "linking_eye.Accordion.ETDRS0", 
         "linking_eye.Accordion.ETDRS4", 
         "linking_eye.Accordion.ETDRS8", 
         "ETDRS_B4", 
         "ETDRS_48", 
         "ETDRS_B8", 
         "Linking.MIND_Accordion.TOTAL_DSC_BLR", 
         "Linking.MIND_Accordion.TOTAL_DSC_F40", 
         "Linking.MIND_Accordion.TOTAL_DSC_F80", 
         "Linking.MIND_Accordion.TOTAL_MMSE_BLR", 
         "Linking.MIND_Accordion.TOTAL_MMSE_F40", 
         "Linking.MIND_Accordion.TOTAL_MMSE_F80", 
         "Linking.MIND_Accordion.STROOP_BLR", 
         "Linking.MIND_Accordion.STROOP_F40", 
         "Linking.MIND_Accordion.STROOP_F80", 
         "Linking.MIND_Accordion.RAVLT_BLR", 
         "Linking.MIND_Accordion.RAVLT_F40", 
         "Linking.MIND_Accordion.RAVLT_F80", 
         "DSST_B40", 
         "DSST_4080", 
         "DSST_B80", 
         "MMSE_B40", 
         "MMSE_4080", 
         "MMSE_B80", 
         "STROOP_B40", 
         "STROOP_4080", 
         "STROOP_B80", 
         "RAVLT_B40", 
         "RAVLT_4080", 
         "RAVLT_B80", 
         "TBV_B40", 
         "WMV_B40", 
         "GMV_B40", 
         "LHC_B40", 
         "RHC_B40", 
         "LLG_B40", 
         "RLG_B40", 
         "TBV_4080", 
         "WMV_4080", 
         "GMV_4080", 
         "LHC_4080", 
         "RHC_4080", 
         "LLG_4080", 
         "RLG_4080", 
         "TBV_B80", 
         "WMV_B80", 
         "GMV_B80", 
         "LHC_B80", 
         "RHC_B80", 
         "LLG_B80", 
         "RLG_B80", 
         "TBV_BLR", 
         "TBV_F40", 
         "TBV_F80", 
         "WM_BLR", 
         "WM_F40", 
         "WM_F80", 
         "GM_BLR", 
         "GM_F40", 
         "GM_F80", 
         "LEFT_HIPPOCAMPUS_BLR", 
         "LEFT_HIPPOCAMPUS_F40", 
         "LEFT_HIPPOCAMPUS_F80", 
         "RIGHT_HIPPOCAMPUS_BLR", 
         "RIGHT_HIPPOCAMPUS_F40", 
         "RIGHT_HIPPOCAMPUS_F80", 
         "Left_LiG___lingual_gyrus_BLR", 
         "Left_LiG___lingual_gyrus_F40", 
         "Left_LiG___lingual_gyrus_F80", 
         "Right_LiG___lingual_gyrus_BLR", 
         "Right_LiG___lingual_gyrus_F40", 
         "Right_LiG___lingual_gyrus_F80",
         "CORPUS_CALLOSUM_BLR",
         "CORPUS_CALLOSUM_F40",
         "CORPUS_CALLOSUM_F80",
         "CC_B40",
         "CC_4080",
         "CC_B80",
         "Left_Ent___entorhinal_area_BLR",
         "Left_Ent___entorhinal_area_F40",
         "Left_Ent___entorhinal_area_F80",
         "Right_Ent___entorhinal_area_BLR",
         "Right_Ent___entorhinal_area_F40",
         "Right_Ent___entorhinal_area_F80",
         "Left_PHG___parahippocampal_gyrus_BLR",
         "Left_PHG___parahippocampal_gyrus_F40",
         "Left_PHG___parahippocampal_gyrus_F80",
         "Right_PHG___parahippocampal_gyrus_BLR",
         "Right_PHG___parahippocampal_gyrus_F40",
         "Right_PHG___parahippocampal_gyrus_F80",
         "Right_ACgG__anterior_cingulate_gyrus_BLR",
         "Right_ACgG__anterior_cingulate_gyrus_F40",
         "Right_ACgG__anterior_cingulate_gyrus_F80",
         "LIMBIC_GM_BLR",
         "LIMBIC_GM_F40",
         "LIMBIC_GM_F80",
         "x2stroke",
         "any_diab_med",
         "TOTAL_BRAIN_VOLUME_ICV_BLR",
         "TOTAL_BRAIN_VOLUME_ICV_F40",
         "TOTAL_BRAIN_VOLUME_ICV_F80",
         "GRAY_MATTER_NORMAL_SUM_BLR",
         "GRAY_MATTER_NORMAL_SUM_F40",
         "GRAY_MATTER_NORMAL_SUM_F80",
         "WHITE_MATTER_NORMAL_SUM_BLR",
         "WHITE_MATTER_NORMAL_SUM_F40",
         "WHITE_MATTER_NORMAL_SUM_F80",
         "WHITE_MATTER_ABNORMAL_SUM_BLR",
         "yjWHITE_MATTER_ABNORMAL_SUM_BLR",
         "thresholdWHITE_MATTER_ABNORMAL_SUM_BLR",
         "WHITE_MATTER_ABNORMAL_SUM_F40",
         "WHITE_MATTER_ABNORMAL_SUM_F80",
         "GRAY_MATTER_ABNORMAL_SUM_BLR",
         "GRAY_MATTER_ABNORMAL_SUM_F40",
         "GRAY_MATTER_ABNORMAL_SUM_F80",
         "TOTAL_BRAIN_VOLUME_ICV_B40",
         "TOTAL_BRAIN_VOLUME_ICV_B80",
         "GRAY_MATTER_NORMAL_SUM_B40",
         "GRAY_MATTER_NORMAL_SUM_B80",
         "WHITE_MATTER_NORMAL_SUM_B40",
         "WHITE_MATTER_NORMAL_SUM_B80",
         "WHITE_MATTER_ABNORMAL_SUM_B40",
         "WHITE_MATTER_ABNORMAL_SUM_B80",
         "GRAY_MATTER_ABNORMAL_SUM_B40",
         "GRAY_MATTER_ABNORMAL_SUM_B80",
  )

retinsadf_limited <- retinsadf[!is.na(retinsadf$linking_eye.Accordion.ETDRS0),]
retinsadf_table1 <- retinsadf_limited
View(retinsadf_limited)

#######################################################################################################

# Table 1: Cohort Characteristics

#######################################################################################################
# Convert Numerically Coded Data to Factors
retinsadf_table1 <- retinsadf_table1 %>%
  mutate(
    sex = factor(sex, levels = c(0, 1), 
                 labels = c("Female", "Male")),
    cvd_hx_baseline = factor(cvd_hx_baseline, levels =c(0, 1), 
                             labels = c("No", "Yes")),
    Linking_concomitantmedsAccordio.MEM_LOSS = factor(Linking_concomitantmedsAccordio.MEM_LOSS, levels =c(0, 1), 
                                                      labels = c("No", "Yes")),
    Linking_concomitantmedsAccordio.ANTI_DEPRESS = factor(Linking_concomitantmedsAccordio.ANTI_DEPRESS, levels =c(0, 1), 
                                                          labels = c("No", "Yes")),
    edu = factor(edu, levels =c(1, 2, 3, 4),
                 labels = c("Less than Highschool", "High School Graduate", "Some College", "College Degree of Higher")),
    cigarett = factor(cigarett, levels =c(1, 2),
                      labels = c("Yes", "No")),
    alcohol = factor(alcohol > 0,
                     levels = c(TRUE, FALSE),
                     labels = c("Yes", "No")),
    x2stroke = factor(x2stroke, levels = c(1, 2), 
                 labels = c("Yes", "No")),
    any_diab_med = factor(any_diab_med, levels = c(0, 1), 
                      labels = c("No", "Yes")),
  )

# Re-Make Non-Log AB Ratio
retinsadf_table1$AB4240_Ratio <- exp(retinsadf_table1$logAb4240_Ratio-1)

# Specify Labels for Table 1
label(retinsadf_table1$baseline_age) <- "Baseline Age (yrs)"
label(retinsadf_table1$race) <- "Race"
label(retinsadf_table1$sex) <- "Sex"
label(retinsadf_table1$edu) <- "Education"
label(retinsadf_table1$yrsdiab) <- "Diabetes Duration (yrs)"
label(retinsadf_table1$linking_hba1cAccordion.HBA1C) <- "Hemoglobin A1c at Baseline (%)"
label(retinsadf_table1$cvd_hx_baseline) <- "CVD History at Baseline"
label(retinsadf_table1$Linking_Lipids.Accordion.CHOL) <- "Cholesterol (mg/dL)"
label(retinsadf_table1$Linking_Lipids.Accordion.LDL) <- "LDL (mg/dL)"
label(retinsadf_table1$Linking_Lipids.Accordion.HDL) <- "HDL (mg/dL)"
label(retinsadf_table1$Linking.otherlabs.Accordion.GFR) <- "GFR (mL/min/1.73 m²)"
label(retinsadf_table1$BMI) <- "BMI (kg/m²)"
label(retinsadf_table1$Linking_blodpressure.Accordion.SBP) <- "Systolic Blood Pressure (mmHg)"
label(retinsadf_table1$Linking_blodpressure.Accordion.DBP) <- "Diastolic Blood Pressure (mmHg)"
label(retinsadf_table1$cigarett) <- "Smoked  in the Last 30 Days"
label(retinsadf_table1$alcohol) <- "Currently Drinks Alcohol"
label(retinsadf_table1$Linking.MIND_Accordion.RAVLT_BLR) <- "RAVLT"
label(retinsadf_table1$Linking.MIND_Accordion.TOTAL_MMSE_BLR) <- "MMSE"
label(retinsadf_table1$Linking.MIND_Accordion.STROOP_BLR) <- "Stroop"
label(retinsadf_table1$Linking.MIND_Accordion.TOTAL_DSC_BLR) <- "DSST"
label(retinsadf_table1$Linking_concomitantmedsAccordio.MEM_LOSS) <- "Any Memory Loss Medication"
label(retinsadf_table1$Linking_concomitantmedsAccordio.ANTI_DEPRESS) <- "Any Anti-Depressant Medication"
label(retinsadf_table1$ptau217) <- "pTau-217 (AU)"
label(retinsadf_table1$AB4240_Ratio) <- "Ab42:40 Ratio"
label(retinsadf_table1$NFL) <- "NFL (pg/mL)"
label(retinsadf_table1$linking_eye.Accordion.ETDRS0) <- "Baseline ETDRS"
label(retinsadf_table1$x2stroke) <- "History of Stroke"
label(retinsadf_table1$any_diab_med) <- "Any Diabetes Medication"

# Stratify Protein Levels
# pTau-217
breaks_ptau <- quantile(retinsadf_table1$logptau217, c(.33, .67, 1))
breaks_ptau <- c(0, breaks_ptau)
labels_ptau <- c('Low pTau-217', 'Medium pTau-217', 'High pTau-217')
levels_ptau <- cut(retinsadf_table1$logptau217, breaks = breaks_ptau, labels = labels_ptau)
retinsadf_table1 <- cbind(retinsadf_table1, levels_ptau)

# Ab42:40 Ratio
breaks_Ab_4240 <- quantile(retinsadf_table1$logAb4240_Ratio, c(.33, .67, 1))
breaks_Ab_4240 <- c(0, breaks_Ab_4240)
labels_Ab_4240 <- c('Low Ab42:40 Ratio', 'Medium Ab42:40 Ratio', 'High Ab42:40 Ratio')
levels_Ab4240 <- cut(retinsadf_table1$logAb4240_Ratio, breaks = breaks_Ab_4240, labels = labels_Ab_4240)
retinsadf_table1 <- cbind(retinsadf_table1, levels_Ab4240)

# NFL
breaks_NFL <- quantile(retinsadf_table1$logNFL, c(.33, .67, 1))
breaks_NFL <- c(0, breaks_NFL)
labels_NFL <- c('Low NFL', 'Medium NFL', 'High NFL')
levels_NFL <- cut(retinsadf_table1$logNFL, breaks = breaks_NFL, labels = labels_NFL)
retinsadf_table1 <- cbind(retinsadf_table1, levels_NFL)

# Generate Table 1
# Specify Statistics to Display
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits = 2), 
       c("",
         "Mean (SD)" = sprintf("%s (&plusmn; %s)", MEAN, SD),
         "Median [IQR]" = sprintf("%s [%s, %s]", MEDIAN, Q1, Q3),
         "Min, Max"  = sprintf("%s, %s", MIN, MAX))
  )
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

table1(~ retinsadf_table1$baseline_age +
         retinsadf_table1$race +
         retinsadf_table1$sex +
         retinsadf_table1$edu +
         retinsadf_table1$yrsdiab +
         retinsadf_table1$linking_hba1cAccordion.HBA1C +
         retinsadf_table1$any_diab_med +
         retinsadf_table1$cvd_hx_baseline +
         retinsadf_table1$x2stroke +
         retinsadf_table1$Linking_Lipids.Accordion.CHOL +
         retinsadf_table1$Linking_Lipids.Accordion.LDL +
         retinsadf_table1$Linking_Lipids.Accordion.HDL +
         retinsadf_table1$Linking.otherlabs.Accordion.GFR +
         retinsadf_table1$BMI +
         retinsadf_table1$Linking_blodpressure.Accordion.SBP +
         retinsadf_table1$Linking_blodpressure.Accordion.DBP +
         retinsadf_table1$cigarett +
         retinsadf_table1$alcohol +
         retinsadf_table1$Linking.MIND_Accordion.RAVLT_BLR +
         retinsadf_table1$Linking.MIND_Accordion.TOTAL_MMSE_BLR +
         retinsadf_table1$Linking.MIND_Accordion.STROOP_BLR +
         retinsadf_table1$Linking.MIND_Accordion.TOTAL_DSC_BLR +
         retinsadf_table1$Linking_concomitantmedsAccordio.MEM_LOSS +
         retinsadf_table1$Linking_concomitantmedsAccordio.ANTI_DEPRESS +
         retinsadf_table1$ptau217 +
         retinsadf_table1$AB4240_Ratio +
         retinsadf_table1$NFL +
         retinsadf_table1$linking_eye.Accordion.ETDRS0, data=retinsadf_table1, render.continuous=my.render.cont, render.categorical=my.render.cat, render.missing=NULL)

table1(~ retinsadf_table1$baseline_age +
         retinsadf_table1$race +
         retinsadf_table1$sex +
         retinsadf_table1$edu +
         retinsadf_table1$yrsdiab +
         retinsadf_table1$linking_hba1cAccordion.HBA1C +
         retinsadf_table1$any_diab_med +
         retinsadf_table1$cvd_hx_baseline +
         retinsadf_table1$x2stroke +
         retinsadf_table1$Linking_Lipids.Accordion.CHOL +
         retinsadf_table1$Linking_Lipids.Accordion.LDL +
         retinsadf_table1$Linking_Lipids.Accordion.HDL +
         retinsadf_table1$Linking.otherlabs.Accordion.GFR +
         retinsadf_table1$BMI +
         retinsadf_table1$Linking_blodpressure.Accordion.SBP +
         retinsadf_table1$Linking_blodpressure.Accordion.DBP +
         retinsadf_table1$cigarett +
         retinsadf_table1$alcohol +
         retinsadf_table1$Linking.MIND_Accordion.TOTAL_MMSE_BLR +
         retinsadf_table1$ptau217 +
         retinsadf_table1$AB4240_Ratio +
         retinsadf_table1$NFL +
         retinsadf_table1$linking_eye.Accordion.ETDRS0 | levels_NFL, data=retinsadf_table1, render.continuous=my.render.cont, render.categorical=my.render.cat, render.missing=NULL)

#######################################################################################################

# Table 2: Cross-Sectional Analysis of Baseline NDD Marker vs ETDRS

#######################################################################################################
# Retinopathy vs Neurodegenerative Markers (Do NDD Markers influence ETDRS scores?)
# ETDRS 0
model1_ETDRS0_ptau <- lm(linking_eye.Accordion.ETDRS0 ~ logptau217, data = retinsadf_limited)
model2_ETDRS0_ptau <- lm(linking_eye.Accordion.ETDRS0~ logptau217 + baseline_age + sex + edu + race, data = retinsadf_limited)
model3_ETDRS0_ptau <- lm(linking_eye.Accordion.ETDRS0 ~ logptau217 + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + retinsadf_limited$cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited)

model1_ETDRS0_NFL <- lm(linking_eye.Accordion.ETDRS0 ~ logNFL, data = retinsadf_limited)
model2_ETDRS0_NFL <- lm(linking_eye.Accordion.ETDRS0~ logNFL + baseline_age + sex + edu + race, data = retinsadf_limited)
model3_ETDRS0_NFL <- lm(linking_eye.Accordion.ETDRS0 ~ logNFL + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited)

model1_ETDRS0_ABRatio <- lm(linking_eye.Accordion.ETDRS0 ~ logAb4240_Ratio, data = retinsadf_limited)
model2_ETDRS0_ABRatio <- lm(linking_eye.Accordion.ETDRS0~ logAb4240_Ratio + baseline_age + sex + edu + race, data = retinsadf_limited)
model3_ETDRS0_ABRatio <- lm(linking_eye.Accordion.ETDRS0 ~ logAb4240_Ratio + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + retinsadf_limited$cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited)

# Model Summaries for Table 2 Retinopathy vs Neurodegenerative Markers
# ETDRS 0
model1_ETDRS0_ptau_list <- list("Model 1"=model1_ETDRS0_ptau,"Model 2"=model2_ETDRS0_ptau,"Model 3"=model3_ETDRS0_ptau)
modelsummary(model1_ETDRS0_ptau_list, statistic = "P-value: {p.value}", stars = TRUE)
model1_ETDRS0_NFL_list <- list("Model 1"=model1_ETDRS0_NFL,"Model 2"=model2_ETDRS0_NFL,"Model 3"=model3_ETDRS0_NFL)
modelsummary(model1_ETDRS0_NFL_list, statistic = "P-value: {p.value}", stars = TRUE)
model1_ETDRS0_ABRatio_list <- list("Model 1"=model1_ETDRS0_ABRatio,"Model 2"=model2_ETDRS0_ABRatio,"Model 3"=model3_ETDRS0_ABRatio)
modelsummary(model1_ETDRS0_ABRatio_list, statistic = "P-value: {p.value}", stars = TRUE)
# ETDRS 0 w/ 95% CI
model1_ETDRS0_ptau_list <- list("Model 1"=model1_ETDRS0_ptau,"Model 2"=model2_ETDRS0_ptau,"Model 3"=model3_ETDRS0_ptau)
modelsummary(model1_ETDRS0_ptau_list, statistic = "conf.int", stars = TRUE)
model1_ETDRS0_NFL_list <- list("Model 1"=model1_ETDRS0_NFL,"Model 2"=model2_ETDRS0_NFL,"Model 3"=model3_ETDRS0_NFL)
modelsummary(model1_ETDRS0_NFL_list, statistic = "conf.int", conf_level = .95, stars = TRUE)
model1_ETDRS0_ABRatio_list <- list("Model 1"=model1_ETDRS0_ABRatio,"Model 2"=model2_ETDRS0_ABRatio,"Model 3"=model3_ETDRS0_ABRatio)
modelsummary(model1_ETDRS0_ABRatio_list, statistic = "conf.int", stars = TRUE)

######################################################################################

# Figure 3: NfL Violin Plot by Clinically Relevant DR Stage

######################################################################################
# Make DR Categories for ETDRS 0
retinsadf_limited <- retinsadf_limited %>%
  mutate(DR_Severity = case_when(
    linking_eye.Accordion.ETDRS0 == 1 ~ "No Retinopathy",
    linking_eye.Accordion.ETDRS0 %in% 2:3 ~ "Mild DR",
    linking_eye.Accordion.ETDRS0 %in% 4:17 ~ "Moderate/Severe DR",
    TRUE ~ NA_character_
  ))

retinsadf_limited <- retinsadf_limited %>%
  mutate(DR_Severity = factor(DR_Severity,
                              levels = c("No Retinopathy", "Mild DR", "Moderate/Severe DR"),
                              ordered = TRUE))

nfl_compare <- compare_means(logNFL ~ DR_Severity,  data = retinsadf_limited,
                             method = "t.test")

nfl_compare_list <- list(
  c("No Retinopathy", "Mild DR"),
  c("No Retinopathy", "Moderate/Severe DR"),
  c("Mild DR", "Moderate/Severe DR"))

nfl_v_plot <- ggplot(data=subset(retinsadf_limited, !is.na(DR_Severity)), aes(x=DR_Severity, y=logNFL, fill = DR_Severity)) + 
  geom_violin()+
  geom_jitter()+
  geom_boxplot(width=0.1)+
  theme_pubclean()+
  labs(x="Diabetic Retinopathy Severity", y = "logNfL (pg/mL)")+
  scale_x_discrete(limits=c("No Retinopathy", "Mild DR", "Moderate/Severe DR"))+
  guides(fill=guide_legend(title="DR Severity:"))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c(
    "No Retinopathy" = "#7CAE00",
    "Mild DR" = "#00BFC4",
    "Moderate/Severe DR" = "#F8766D"
  ))+
  stat_compare_means(method = "anova", label.y = 1.7, label.x = 2.75)+
  stat_compare_means(comparisons = nfl_compare_list, label = "p.format", method = "t.test",
                     tip.length = 0, step.increase = 0.06, label.y = 4.2,
                     ref.group = "No Retinopathy")
nfl_v_plot

#######################################################################################################

# Table 4: Forest Plot

#######################################################################################################
# Re-level Variables
retinsadf_table1$sex <- relevel(retinsadf_table1$sex, ref = "Female")
retinsadf_table1$alcohol <- relevel(retinsadf_table1$alcohol, ref = "No")
retinsadf_table1$x2stroke <- relevel(retinsadf_table1$x2stroke, ref = "No")
retinsadf_table1$cvd_hx_baseline <- relevel(retinsadf_table1$cvd_hx_baseline, ref = "No")
retinsadf_table1$cigarett <- relevel(retinsadf_table1$cigarett, ref = "No")
retinsadf_table1$any_diab_med <- relevel(retinsadf_table1$any_diab_med, ref = "No")

# Create List of Categorical Effect Modifiers
cat_effect_modifiers <- c("sex",
                          "edu",
                          "race",
                          "any_diab_med",
                          "cigarett",
                          "cvd_hx_baseline",
                          "alcohol",
                          "x2stroke"
)

# Make Categorical Effect Modifiers Factors
retinsadf_table1 <- retinsadf_table1 %>%
  mutate(
    across(
      all_of(cat_effect_modifiers),
      as.factor
    )
  )
View(retinsadf_table1)

cat_forest_results <- map_dfr(
  cat_effect_modifiers,
  function(modifier) {
    
    formula <- as.formula(
      paste(
        "linking_eye.Accordion.ETDRS0 ~ logNFL *",
        modifier
      )
    )
    
    model <- lm(
      formula,
      data = retinsadf_table1
    )
    
    tidy(
      model,
      conf.int = TRUE
    ) %>%
      filter(
        grepl("^logNFL:", term)
      ) %>%
      mutate(
        Modifier = modifier
      )
  }
)
cat_forest_results

# Make Labels
cat_forest_results <- cat_forest_results %>%
  mutate(
    label = case_when(
      term == "logNFL:sexMale" ~ "Male vs Female",
      grepl("logNFL:eduHigh School Graduate", term) ~
        "High School Graduate vs < High School",
      grepl("logNFL:eduSome College", term) ~
        "Some College vs < High School",
      grepl("logNFL:eduCollege Degree of Higher", term) ~
        "College Degree or Higher vs < High School",
      grepl("logNFL:raceBlack", term) ~
        "Black vs White",
      grepl("logNFL:raceHispanic", term) ~
        "Hispanic vs White",
      grepl("logNFL:raceOther", term) ~
        "Other vs White",
      grepl("logNFL:any_diab_medYes", term) ~
        "Diabetes medication: Yes vs No",
      grepl("logNFL:cigarettYes", term) ~
        "Smoking: Yes vs No",
      grepl("logNFL:cvd_hx_baselineYes", term) ~
        "CVD history: Yes vs No",
      grepl("logNFL:alcoholYes", term) ~
        "Alcohol: Yes vs No",
      grepl("logNFL:x2strokeYes", term) ~
        "Stroke: Yes vs No",
      TRUE ~ term
    )
  )
cat_forest_results

# Data Formatting
cat_plot_data <- cat_forest_results %>%
  mutate(
    term_factor = factor(label, levels = rev(unique(label))),
    p_formatted = ifelse(p.value < 0.05, "<0.05", sprintf("%.3f", p.value)),
    est_ci_text = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)
  )

# Forest Plot 
cat_forest_plot <- ggplot(cat_plot_data, aes(y = term_factor, x = estimate)) +
  geom_point(shape = 18, size = 3) +  
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_discrete(name = "") +
  labs(
    x = "Beta (95% CI)",
    y = NULL,
    title = "Categorical Effect Modifiers of logNfL ~ ETDRS DRSS"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.placement = "outside", # Places modifier titles to the far left
    strip.text.y.left = element_text(angle = 0, face = "bold"), # Keep headers horizontal
    panel.spacing = unit(0.5, "lines"),
    axis.text.y = element_text(size = 10)
  )
cat_forest_plot

# Data Table
cat_fp_table <- ggplot(cat_plot_data, aes(y = term_factor)) +
  # Beta & 95% CI column
  geom_text(aes(x = 1, label = est_ci_text), hjust = 0.5, size = 3.5) +
  # P-value column
  geom_text(aes(x = 2, label = p_formatted), hjust = 0.5, size = 3.5) +
  # Column headers
  annotate("text", x = 1, y = length(unique(cat_plot_data$label)) + 0.8, 
           label = "Beta (95% CI)", fontface = "bold", size = 3.8) +
  annotate("text", x = 2, y = length(unique(cat_plot_data$label)) + 0.8, 
           label = "P-value", fontface = "bold", size = 3.8) +
  scale_x_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) +
  # Match y-axis bounds precisely so rows align
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(t = 30, r = 10, b = 30, l = 10)
  )

# Plot Forest Plot & Table Side by Side
cat_combined_plot <- cat_forest_plot + cat_fp_table + plot_layout(widths = c(2.5, 1))

# Display plot
cat_combined_plot

# Save
ggsave(
  filename = "cat_nfl_forest_plot_combined.pdf",
  plot = cat_combined_plot,
  width = 10,       # Width in inches
  height = 6,       # Adjust height based on number of modifier rows
  units = "in",
  dpi = 300
)


# Make List on Continuous Effect Modifiers
cont_effect_modifiers <- c("baseline_age",
                           "yrsdiab",
                           "linking_hba1cAccordion.HBA1C",
                           "Linking_blodpressure.Accordion.SBP",
                           "BMI",
                           "Linking_Lipids.Accordion.CHOL",
                           "Linking.otherlabs.Accordion.GFR",
                           "Linking_Lipids.Accordion.LDL",
                           "Linking_Lipids.Accordion.HDL",
                           "Linking.MIND_Accordion.TOTAL_MMSE_BLR",
                           "TOTAL_BRAIN_VOLUME_ICV_BLR",
                           "GRAY_MATTER_NORMAL_SUM_BLR",
                           "WHITE_MATTER_NORMAL_SUM_BLR",
                           "logWHITE_MATTER_ABNORMAL_SUM_BLR"
)

cont_forest_labels <- c(
  baseline_age = "Age (yrs)",
  yrsdiab = "Diabetes Duration (yrs)",
  linking_hba1cAccordion.HBA1C = "Hemoglobin A1c (%)",
  Linking_blodpressure.Accordion.SBP = "Systolic BP (mmHg)",
  BMI = "BMI (kg/m²)",
  Linking_Lipids.Accordion.CHOL = "Total cholesterol (mg/dL)",
  Linking.otherlabs.Accordion.GFR = "eGFR (mL/min/1.73 m²)",
  Linking_Lipids.Accordion.LDL = "LDL (mg/dL)",
  Linking_Lipids.Accordion.HDL = "HDL (mg/dL)",
  Linking.MIND_Accordion.TOTAL_MMSE_BLR = "MMSE",
  TOTAL_BRAIN_VOLUME_ICV_BLR = "Total Brain Volume",
  GRAY_MATTER_NORMAL_SUM_BLR = "Gray Matter Volume",
  WHITE_MATTER_NORMAL_SUM_BLR = "White Matter Volume",
  logWHITE_MATTER_ABNORMAL_SUM_BLR = "Abnormal White Matter Volume"
)


cont_forest_results <- map_dfr(
  cont_effect_modifiers,
  function(modifier) {
    
    formula <- as.formula(
      paste(
        "linking_eye.Accordion.ETDRS0 ~ logNFL *",
        modifier
      )
    )
    
    model <- lm(
      formula,
      data = retinsadf_table1
    )
    
    tidy(
      model,
      conf.int = TRUE
    ) %>%
      filter(
        grepl("^logNFL:", term)
      ) %>%
      mutate(
        Modifier = modifier
      )
  }
)
cont_forest_results
# Make Labels
retinsadf_table1$linking_hba1cAccordion.HBA1C
retinsadf_table1$yrsdiab
cont_forest_results <- cont_forest_results %>%
  mutate(
    label = case_when(
      term == "logNFL:baseline_age" ~ "Age (yrs)",
      grepl("logNFL:linking_hba1cAccordion.HBA1C", term) ~
        "Hemoglobin A1c (%)",
      grepl("logNFL:yrsdiab", term) ~
        "Diabetes Duration (yrs)",
      grepl("logNFL:Linking_blodpressure.Accordion.SBP", term) ~
        "Systolic BP (mmHg)",
      grepl("logNFL:BMI", term) ~
        "BMI (kg/m²)",
      grepl("logNFL:Linking_Lipids.Accordion.CHOL", term) ~
        "Total cholesterol (mg/dL)",
      grepl("logNFL:Linking.otherlabs.Accordion.GFR", term) ~
        "GFR (mL/min/1.73 m²)",
      grepl("logNFL:Linking_Lipids.Accordion.LDL", term) ~
        "LDL (mg/dL)",
      grepl("logNFL:Linking_Lipids.Accordion.HDL", term) ~
        "HDL (mg/dL)",
      grepl("logNFL:Linking.MIND_Accordion.TOTAL_MMSE_BLR", term) ~
        "MMSE",
      grepl("logNFL:TOTAL_BRAIN_VOLUME_ICV_BLR", term) ~
        "Total Brain Volume",
      grepl("logNFL:GRAY_MATTER_NORMAL_SUM_BLR", term) ~
        "Gray Matter Volume",
      grepl("logNFL:WHITE_MATTER_NORMAL_SUM_BLR", term) ~
        "White Matter Volume",
      grepl("logNFL:logWHITE_MATTER_ABNORMAL_SUM_BLR", term) ~
        "Abnormal White Matter Volume", TRUE ~ term
    )
  )
cont_forest_results

# Data Formatting
cont_plot_data <- cont_forest_results %>%
  mutate(
    term_factor = factor(label, levels = rev(unique(label))),
    p_formatted = ifelse(p.value < 0.05, "<0.05", sprintf("%.3f", p.value)),
    est_ci_text = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)
  )

# Forest Plot 
cont_forest_plot <- ggplot(cont_plot_data, aes(y = term_factor, x = estimate)) +
  geom_point(shape = 18, size = 3) +  
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_discrete(name = "") +
  labs(
    x = "Beta (95% CI)",
    y = NULL,
    title = "Continuous Effect Modifiers of logNfL ~ ETDRS DRSS"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.placement = "outside", # Places modifier titles to the far left
    strip.text.y.left = element_text(angle = 0, face = "bold"), # Keep headers horizontal
    panel.spacing = unit(0.5, "lines"),
    axis.text.y = element_text(size = 10)
  )
cont_forest_plot

# Data Table
cont_fp_table <- ggplot(cont_plot_data, aes(y = term_factor)) +
  # Beta & 95% CI column
  geom_text(aes(x = 1, label = est_ci_text), hjust = 0.5, size = 3.5) +
  # P-value column
  geom_text(aes(x = 2, label = p_formatted), hjust = 0.5, size = 3.5) +
  # Column headers
  annotate("text", x = 1, y = length(unique(cont_plot_data$label)) + 0.8, 
           label = "Beta (95% CI)", fontface = "bold", size = 3.8) +
  annotate("text", x = 2, y = length(unique(cont_plot_data$label)) + 0.8, 
           label = "P-value", fontface = "bold", size = 3.8) +
  scale_x_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) +
  # Match y-axis bounds precisely so rows align
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(t = 30, r = 10, b = 30, l = 10)
  )

#retinsadf_table1$GRAY_MATTER_NORMAL_SUM_BLR
#test_model <- lm(linking_eye.Accordion.ETDRS0 ~ logNFL * GRAY_MATTER_NORMAL_SUM_BLR, data = retinsadf_table1)
#modelsummary(test_model, statistic = "P-value: {p.value}", stars = TRUE)

# Plot Forest Plot & Table Side by Side
cont_combined_plot <- cont_forest_plot + cont_fp_table + plot_layout(widths = c(2.5, 1))

# Display plot
cont_combined_plot

# Save
ggsave(
  filename = "cont_nfl_forest_plot_combined.pdf",
  plot = cont_combined_plot,
  width = 10,       # Width in inches
  height = 6,       # Adjust height based on number of modifier rows
  units = "in",
  dpi = 300
)

#######################################################################################################

# Additional Analyses to Check Robustness of Findings

#######################################################################################################
# Created Ordered Category for ETDRS

# ETDRS Scale
# 1: None
# 2-3: Mild
# 4-17: Moderate/Severe
retinsadf_limited <- retinsadf_limited %>%
  mutate(DR_Severity_OLR = case_when(
      linking_eye.Accordion.ETDRS0 == 1 ~ "No_DR",
      linking_eye.Accordion.ETDRS0 %in% 2:3 ~ "Mild_DR",
      linking_eye.Accordion.ETDRS0 %in% 4:17 ~ "Mod_Severe_DR",
    TRUE ~ NA_character_
  ))

# Order DR Severity
retinsadf_limited$DR_Severity_OLR <- factor(
  retinsadf_limited$DR_Severity_OLR,
  levels = c("No_DR", "Mild_DR", "Mod_Severe_DR"),
  ordered = TRUE
)

# Check that severity is ordered and what the order is
is.factor(retinsadf_limited$DR_Severity_OLR)
is.ordered(retinsadf_limited$DR_Severity_OLR)
levels(retinsadf_limited$DR_Severity_OLR)

# Count the number of observations in each category
retinsadf_limited |> count(DR_Severity_OLR)

# Ordinal Logistic Regression for Categorical Outcomes
set.seed(123)
remove.packages("MASS")
install.packages("MASS")
library(MASS) # If MASS fails to load, uninstall and restart R session
polr_model1_ETDRS0_ptau <- polr(DR_Severity_OLR ~ logptau217, data = retinsadf_limited, Hess = TRUE)
polr_model2_ETDRS0_ptau <- polr(DR_Severity_OLR ~ logptau217 + baseline_age + sex + edu + race, data = retinsadf_limited, Hess = TRUE)
polr_model3_ETDRS0_ptau <- polr(DR_Severity_OLR ~ logptau217 + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited, Hess = TRUE)

polr_model1_ETDRS0_NFL <- polr(DR_Severity_OLR ~ logNFL, data = retinsadf_limited, Hess = TRUE)
polr_model2_ETDRS0_NFL <- polr(DR_Severity_OLR ~ logNFL + baseline_age + sex + edu + race, data = retinsadf_limited, Hess = TRUE)
polr_model3_ETDRS0_NFL <- polr(DR_Severity_OLR ~ logNFL + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited, Hess = TRUE)

polr_model1_ETDRS0_ABRatio <- polr(DR_Severity_OLR ~ logAb4240_Ratio, data = retinsadf_limited, Hess = TRUE)
polr_model2_ETDRS0_ABRatio <- polr(DR_Severity_OLR ~ logAb4240_Ratio + baseline_age + sex + edu + race, data = retinsadf_limited, Hess = TRUE)
polr_model3_ETDRS0_ABRatio <- polr(DR_Severity_OLR ~ logAb4240_Ratio + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited, Hess = TRUE)

# Model Summaries for Ordinal Logistic Regression of Retinopathy vs Neurodegenerative Markers
polr_model1_ETDRS0_ptau_list <- list("Model 1"=polr_model1_ETDRS0_ptau,"Model 2"=model2_ETDRS0_ptau,"Model 3"=model3_ETDRS0_ptau)
modelsummary(polr_model1_ETDRS0_ptau_list, statistic = "P-value: {p.value}", stars = TRUE)
polr_model1_ETDRS0_NFL_list <- list("Model 1"=polr_model1_ETDRS0_NFL,"Model 2"=model2_ETDRS0_NFL,"Model 3"=model3_ETDRS0_NFL)
modelsummary(polr_model1_ETDRS0_NFL_list, statistic = "P-value: {p.value}", stars = TRUE)
polr_model1_ETDRS0_ABRatio_list <- list("Model 1"=polr_model1_ETDRS0_ABRatio,"Model 2"=model2_ETDRS0_ABRatio,"Model 3"=model3_ETDRS0_ABRatio)
modelsummary(polr_model1_ETDRS0_ABRatio_list, statistic = "P-value: {p.value}", stars = TRUE)
# ETDRS 0 w/ 95% CI
polr_model1_ETDRS0_ptau_list <- list("Model 1"=polr_model1_ETDRS0_ptau,"Model 2"=model2_ETDRS0_ptau,"Model 3"=model3_ETDRS0_ptau)
modelsummary(polr_model1_ETDRS0_ptau_list, statistic = "conf.int", stars = TRUE)
polr_model1_ETDRS0_NFL_list <- list("Model 1"=polr_model1_ETDRS0_NFL,"Model 2"=model2_ETDRS0_NFL,"Model 3"=model3_ETDRS0_NFL)
modelsummary(polr_model1_ETDRS0_NFL_list, statistic = "conf.int", conf_level = .95, stars = TRUE)
polr_model1_ETDRS0_ABRatio_list <- list("Model 1"=polr_model1_ETDRS0_ABRatio,"Model 2"=model2_ETDRS0_ABRatio,"Model 3"=model3_ETDRS0_ABRatio)
modelsummary(polr_model1_ETDRS0_ABRatio_list, statistic = "conf.int", stars = TRUE)

# Check Variance Inflation Factor; currently set for model 3
library(car)
lm_check_PT <- lm(as.numeric(DR_Severity_OLR) ~ logptau217  + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited)
vif(lm_check_PT)
lm_check_NFL <- lm(as.numeric(DR_Severity_OLR) ~ logNFL + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited)
vif(lm_check_NFL)
lm_check_AB <- lm(as.numeric(DR_Severity_OLR) ~ logAb4240_Ratio + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + cvd_hx_baseline + Linking.otherlabs.Accordion.GFR, data = retinsadf_limited)
vif(lm_check_AB)
retinsadf_limited |> count(cvd_hx_baseline)

# Check Proportional Odds of OLR using Brant Test
install.packages("brant")
library(brant)
brant(polr_model1_ETDRS0_ptau)
brant(polr_model2_ETDRS0_ptau)
brant(polr_model3_ETDRS0_ptau)

brant(polr_model1_ETDRS0_NFL)
brant(polr_model2_ETDRS0_NFL)
brant(polr_model3_ETDRS0_NFL)

brant(polr_model1_ETDRS0_ABRatio)
brant(polr_model2_ETDRS0_ABRatio)
brant(polr_model3_ETDRS0_ABRatio)

# Confidence Intervals
cip1 <- confint.default(polr_model1_ETDRS0_ptau)
cip2 <- confint.default(polr_model2_ETDRS0_ptau)
cip3 <- confint.default(polr_model3_ETDRS0_ptau)

cin1 <- confint.default(polr_model1_ETDRS0_NFL)
cin2 <- confint.default(polr_model2_ETDRS0_NFL)
cin3 <- confint.default(polr_model3_ETDRS0_NFL)

cia1 <- confint.default(polr_model1_ETDRS0_ABRatio)
cia2 <- confint.default(polr_model2_ETDRS0_ABRatio)
cia3 <- confint.default(polr_model3_ETDRS0_ABRatio)

# Odds Ratios
exp(coef(polr_model1_ETDRS0_ptau))
exp(coef(polr_model2_ETDRS0_ptau))
exp(coef(polr_model3_ETDRS0_ptau))

exp(coef(polr_model1_ETDRS0_NFL))
exp(coef(polr_model2_ETDRS0_NFL))
exp(coef(polr_model3_ETDRS0_NFL))

exp(coef(polr_model1_ETDRS0_ABRatio))
exp(coef(polr_model2_ETDRS0_ABRatio))
exp(coef(polr_model3_ETDRS0_ABRatio))

# Output CIs + ORs
exp(cbind(OR = coef(polr_model1_ETDRS0_ptau), cip1))
exp(cbind(OR = coef(polr_model2_ETDRS0_ptau), cip2))
exp(cbind(OR = coef(polr_model3_ETDRS0_ptau), cip3))

exp(cbind(OR = coef(polr_model1_ETDRS0_NFL), cin1))
exp(cbind(OR = coef(polr_model2_ETDRS0_NFL), cin2))
exp(cbind(OR = coef(polr_model3_ETDRS0_NFL), cin3))

exp(cbind(OR = coef(polr_model1_ETDRS0_ABRatio), cia1))
exp(cbind(OR = coef(polr_model2_ETDRS0_ABRatio), cia2))
exp(cbind(OR = coef(polr_model3_ETDRS0_ABRatio), cia3))

# NfL does not meet proportional odds assumption for ordinal linear regression

# Perform partial proportional odds model that allows NfL to vary across levels
# without forcing it to be proportional. Only running for NfL as others did not
# violate proportional odds assumption

# PPO using vglm from VGAM package
install.packages("VGAM")
library(VGAM)

ppo_model1_ETDRS0_NFL <- vglm(
  DR_Severity_OLR ~ logNFL,
  family = cumulative(parallel = FALSE ~ logNFL), 
  data = retinsadf_limited
)
summary(ppo_model1_ETDRS0_NFL)

ppo_model2_ETDRS0_NFL <- vglm(
  DR_Severity_OLR ~ logNFL + baseline_age + sex + edu + race,
  family = cumulative(parallel = FALSE ~ logNFL), 
  data = retinsadf_limited
)
summary(ppo_model2_ETDRS0_NFL)

ppo_model3_ETDRS0_NFL <- vglm(
  DR_Severity_OLR ~ logNFL + baseline_age + sex + edu + race + cigarett + Linking_blodpressure.Accordion.SBP + BMI + Linking_Lipids.Accordion.CHOL + cvd_hx_baseline + Linking.otherlabs.Accordion.GFR,
  family = cumulative(parallel = FALSE ~ logNFL), 
  data = retinsadf_limited
)
summary(ppo_model3_ETDRS0_NFL)





