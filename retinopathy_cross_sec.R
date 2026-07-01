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
# Subset Data for Forestploter
fp_retin <- retinsadf_limited %>%
  mutate(
    sex = factor(sex, levels = c(0, 1), 
                 labels = c("Female", "Male")),
    cvd_hx_baseline = factor(cvd_hx_baseline, levels =c(0, 1), 
                             labels = c("No", "Yes")),
    edu = factor(
      ifelse(edu %in% c(1, 2),
             "No College",
             "Some College or Above"),
      levels = c("No College", "Some College or Above")),
    cigarett = factor(cigarett, levels =c(1, 2),
                      labels = c("Yes", "No")),
    alcohol = factor(alcohol > 0,
                     levels = c(TRUE, FALSE),
                     labels = c("Yes", "No")),
    x2stroke = factor(x2stroke, levels = c(1, 2), 
                      labels = c("Yes", "No")),
    any_diab_med = factor(any_diab_med, levels = c(0, 1), 
                          labels = c("No", "Yes")),
    age_group = ifelse(
      baseline_age < median(baseline_age, na.rm = TRUE),
      paste0("<", round(median(baseline_age, na.rm = TRUE), 1)),
      paste0("≥", round(median(baseline_age, na.rm = TRUE), 1))),
diabetes_duration = ifelse(
      yrsdiab < median(yrsdiab, na.rm = TRUE),
      paste0("<", round(median(yrsdiab, na.rm = TRUE), 1)),
      paste0("≥", round(median(yrsdiab, na.rm = TRUE), 1))),
a1c_group = ifelse(
      linking_hba1cAccordion.HBA1C < median(linking_hba1cAccordion.HBA1C, na.rm = TRUE),
      paste0("<", round(median(linking_hba1cAccordion.HBA1C, na.rm = TRUE), 1)),
      paste0("≥", round(median(linking_hba1cAccordion.HBA1C, na.rm = TRUE), 1))),
cholesterol = ifelse(
      Linking_Lipids.Accordion.CHOL < median(Linking_Lipids.Accordion.CHOL, na.rm = TRUE),
      paste0("<", round(median(Linking_Lipids.Accordion.CHOL, na.rm = TRUE), 1)),
      paste0("≥", round(median(Linking_Lipids.Accordion.CHOL, na.rm = TRUE), 1))),
LDL = ifelse(
      Linking_Lipids.Accordion.LDL < median(Linking_Lipids.Accordion.LDL, na.rm = TRUE),
      paste0("<", round(median(Linking_Lipids.Accordion.LDL, na.rm = TRUE), 1)),
      paste0("≥", round(median(Linking_Lipids.Accordion.LDL, na.rm = TRUE), 1))),
HDL = ifelse(
      Linking_Lipids.Accordion.HDL < median(Linking_Lipids.Accordion.HDL, na.rm = TRUE),
      paste0("<", round(median(Linking_Lipids.Accordion.HDL, na.rm = TRUE), 1)),
      paste0("≥", round(median(Linking_Lipids.Accordion.HDL, na.rm = TRUE), 1))),
GFR = ifelse(
      Linking.otherlabs.Accordion.GFR < median(Linking.otherlabs.Accordion.GFR, na.rm = TRUE),
      paste0("<", round(median(Linking.otherlabs.Accordion.GFR, na.rm = TRUE), 1)),
      paste0("≥", round(median(Linking.otherlabs.Accordion.GFR, na.rm = TRUE), 1))),
BMI = ifelse(
     BMI < median(BMI, na.rm = TRUE),
      paste0("<", round(median(BMI, na.rm = TRUE), 1)),
      paste0("≥", round(median(BMI, na.rm = TRUE), 1))),
SBP = ifelse(
      Linking_blodpressure.Accordion.SBP < median(Linking_blodpressure.Accordion.SBP, na.rm = TRUE),
      paste0("<", round(median(Linking_blodpressure.Accordion.SBP, na.rm = TRUE), 1)),
      paste0("≥", round(median(Linking_blodpressure.Accordion.SBP, na.rm = TRUE), 1))),
DBP = ifelse(
      Linking_blodpressure.Accordion.DBP < median(Linking_blodpressure.Accordion.DBP, na.rm = TRUE),
      paste0("<", round(median(Linking_blodpressure.Accordion.DBP, na.rm = TRUE), 1)),
      paste0("≥", round(median(Linking_blodpressure.Accordion.DBP, na.rm = TRUE), 1))),
MMSE = ifelse(
      Linking.MIND_Accordion.TOTAL_MMSE_BLR < median(Linking.MIND_Accordion.TOTAL_MMSE_BLR, na.rm = TRUE),
      paste0("<", round(median(Linking.MIND_Accordion.TOTAL_MMSE_BLR, na.rm = TRUE), 1)),
      paste0("≥", round(median(Linking.MIND_Accordion.TOTAL_MMSE_BLR, na.rm = TRUE), 1))),
TBV = ifelse(
      TOTAL_BRAIN_VOLUME_ICV_BLR < median(TOTAL_BRAIN_VOLUME_ICV_BLR, na.rm = TRUE),
      paste0("<", round(median(TOTAL_BRAIN_VOLUME_ICV_BLR, na.rm = TRUE), 1)),
      paste0("≥", round(median(TOTAL_BRAIN_VOLUME_ICV_BLR, na.rm = TRUE), 1))),
GMV = ifelse(
      GRAY_MATTER_NORMAL_SUM_BLR < median(GRAY_MATTER_NORMAL_SUM_BLR, na.rm = TRUE),
      paste0("<", round(median(GRAY_MATTER_NORMAL_SUM_BLR, na.rm = TRUE), 1)),
      paste0("≥", round(median(GRAY_MATTER_NORMAL_SUM_BLR, na.rm = TRUE), 1))),
WMV = ifelse(
      WHITE_MATTER_NORMAL_SUM_BLR < median(WHITE_MATTER_NORMAL_SUM_BLR, na.rm = TRUE),
      paste0("<", round(median(WHITE_MATTER_NORMAL_SUM_BLR, na.rm = TRUE), 1)),
      paste0("≥", round(median(WHITE_MATTER_NORMAL_SUM_BLR, na.rm = TRUE), 1))),
WMLV = ifelse(
      WHITE_MATTER_ABNORMAL_SUM_BLR < median(WHITE_MATTER_ABNORMAL_SUM_BLR, na.rm = TRUE),
      paste0("<", round(median(WHITE_MATTER_ABNORMAL_SUM_BLR, na.rm = TRUE), 1)),
      paste0("≥", round(median(WHITE_MATTER_ABNORMAL_SUM_BLR, na.rm = TRUE), 1))),
  )

run_subgroup <- function(data, subgroup_var, subgroup_name) {
  
  if (!subgroup_var %in% names(data)) {
    stop(paste("Variable not found:", subgroup_var))
  }
  
  levels_sub <- na.omit(unique(data[[subgroup_var]]))
  
  results <- lapply(levels_sub, function(level) {
    
    subset_data <- data %>%
      filter(.data[[subgroup_var]] == level) %>%
      filter(
        !is.na(linking_eye.Accordion.ETDRS0),
        !is.na(logNFL)
      )
    
    if (nrow(subset_data) < 2) {
      return(NULL)
    }
    
    model <- lm(linking_eye.Accordion.ETDRS0 ~ logNFL,
                data = subset_data)
    
    broom::tidy(model) %>%
      filter(term == "logNFL") %>%
      mutate(
        subgroup = subgroup_name,
        level = as.character(level),
        n_obs = stats::nobs(model)
      )
  })
  
  bind_rows(results)
}

interaction_p <- function(data, subgroup_var) {
  
  model_data <- data %>%
    filter(
      !is.na(.data[[subgroup_var]]),
      !is.na(linking_eye.Accordion.ETDRS0),
      !is.na(logNFL)
    )
  
  if (nrow(model_data) < 3) {
    return(NA_real_)
  }
  
  formula_txt <- paste(
    "linking_eye.Accordion.ETDRS0 ~ logNFL *",
    subgroup_var
  )
  
  model <- lm(as.formula(formula_txt), data = model_data)
  
  pvals <- coef(summary(model))[
    grep("logNFL:", rownames(coef(summary(model)))),
    "Pr(>|t|)"
  ]
  
  if (length(pvals) == 0) NA_real_ else pvals[1]
}

forest_data <- bind_rows(
  run_subgroup(fp_retin, "age_group", "Age"),
  run_subgroup(fp_retin, "sex", "Sex"),
  run_subgroup(fp_retin, "edu", "Education"),
  run_subgroup(fp_retin, "diabetes_duration", "Diabetes Duration"),
  run_subgroup(fp_retin, "a1c_group", "HbA1c"),
  run_subgroup(fp_retin, "any_diab_med", "Any T2DM Medication"),
  run_subgroup(fp_retin, "cvd_hx_baseline", "CVD History"),
  run_subgroup(fp_retin, "x2stroke", "Stroke History"),
  run_subgroup(fp_retin, "cholesterol", "Total Cholesterol"),
  run_subgroup(fp_retin, "LDL", "LDL"),
  run_subgroup(fp_retin, "HDL", "HDL"),
  run_subgroup(fp_retin, "GFR", "GFR"),
  run_subgroup(fp_retin, "BMI", "BMI"),
  run_subgroup(fp_retin, "cigarett", "Currently Smokes"),
  run_subgroup(fp_retin, "alcohol", "Currently Drinks"),
  run_subgroup(fp_retin, "MMSE", "MMSE"),
  run_subgroup(fp_retin, "TBV", "TBV"),
  run_subgroup(fp_retin, "WMV", "WMV"),
  run_subgroup(fp_retin, "GMV", "GMV"),
  run_subgroup(fp_retin, "WMLV", "WMLV")
)

interaction_tbl <- tibble(
  subgroup = c("Age",
               "Sex",
               "Education",
               "Diabetes Duration",
               "HbA1c",
               "Any T2DM Medication",
               "CVD History",
               "Stroke History",
               "Total Cholesterol",
               "LDL",
               "HDL",
               "GFR",
               "BMI",
               "Currently Smokes",
               "Currently Drinks",
               "MMSE",
               "TBV",
               "WMV",
               "GMV",
               "WMLV"
               ),
p_interaction = c(
    interaction_p(fp_retin, "age_group"),
    interaction_p(fp_retin, "sex"),
    interaction_p(fp_retin, "edu"),
    interaction_p(fp_retin, "diabetes_duration"),
    interaction_p(fp_retin, "a1c_group"),
    interaction_p(fp_retin, "any_diab_med"),
    interaction_p(fp_retin, "cvd_hx_baseline"),
    interaction_p(fp_retin, "x2stroke"),
    interaction_p(fp_retin, "cholesterol"),
    interaction_p(fp_retin, "LDL"),
    interaction_p(fp_retin, "HDL"),
    interaction_p(fp_retin, "GFR"),
    interaction_p(fp_retin, "BMI"),
    interaction_p(fp_retin, "cigarett"),
    interaction_p(fp_retin, "alcohol"),
    interaction_p(fp_retin, "MMSE"),
    interaction_p(fp_retin, "TBV"),
    interaction_p(fp_retin, "WMV"),
    interaction_p(fp_retin, "GMV"),
    interaction_p(fp_retin, "WMLV")
  )
)

forest_data <- forest_data %>%
  left_join(interaction_tbl, by = "subgroup") %>%
  mutate(
    subgroup = factor(
      subgroup,
      levels = c(
        "Age",
        "Sex",
        "Education",
        "Diabetes Duration",
        "HbA1c",
        "Any T2DM Medication",
        "CVD History",
        "Stroke History",
        "Total Cholesterol",
        "LDL",
        "HDL",
        "GFR",
        "BMI",
        "Currently Smokes",
        "Currently Drinks",
        "MMSE",
        "TBV",
        "WMV",
        "GMV",
        "WMLV"
      )
    )
  ) %>%
  arrange(subgroup, level) %>%
  group_by(subgroup) %>%
  mutate(
    display_subgroup = if_else(
      row_number() == 1,
      as.character(subgroup),
      ""
    ),
    p_int_display = if_else(
      row_number() == 1,
      if_else(
        is.na(p_interaction),
        "",
        if_else(
          p_interaction < 0.001,
          "<0.001",
          sprintf("%.3f", p_interaction)
        )
      ),
      ""
    )
  ) %>%
  ungroup() %>%
  mutate(
    display_level = as.character(level),
    beta_ci = sprintf(
      "%.2f (%.2f, %.2f)",
      estimate,
      estimate - 1.96 * std.error,
      estimate + 1.96 * std.error
    ),
    p_value_clean = if_else(
      p.value < 0.001,
      "<0.001",
      sprintf("%.3f", p.value)
    )
  )

tabletext <- cbind(
  c("Subgroup", forest_data$display_subgroup),
  c("Level", forest_data$display_level),
  c("N", forest_data$n_obs),
  c("β (95% CI)", forest_data$beta_ci),
  c("P value", forest_data$p_value_clean),
  c("P interaction", forest_data$p_int_display)
)

forestplot(
  labeltext = tabletext,
  mean  = c(NA, forest_data$estimate),
  lower = c(NA, forest_data$estimate - 1.96 * forest_data$std.error),
  upper = c(NA, forest_data$estimate + 1.96 * forest_data$std.error),
  zero = 0,
  boxsize = 0.20,
  lineheight = unit(8, "mm"),
  col = fpColors(
    box = "black",
    line = "black",
    zero = "gray50"
  ),
  xlab = "Beta Coefficient (95% CI)",
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 12),
    summary = gpar(fontface = "bold", fontsize = 12), # Bolds the Subgroup headers
    ticks = gpar(fontsize = 12),
    xlab  = gpar(fontsize = 16),
  ),
  hrzl_lines = list(
    "1" = gpar(lwd = 1.5),
    "2" = gpar(lwd = 1)
  ),
  new_page = TRUE
)
hist(retinsadf_limited$linking_eye.Accordion.ETDRS0)
retinsadf_limited |> count(linking_eye.Accordion.ETDRS8)

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





