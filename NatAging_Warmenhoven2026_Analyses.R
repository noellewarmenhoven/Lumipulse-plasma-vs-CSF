

#######Codes for analyses in Warmenhoven et al. (2026) Nat Aging#######
# If using, please give appropriate credits. 


#Load libraries----
library(readxl)
library(tidyverse)
library(pROC)
library(boot)
library(gtsummary)
library(ggpubr)
library(ggbeeswarm)
library(pROC)
library(see)

#Load datasets----

##BF
bf <- read_xlsx("~/Documents/Projects/LP vs CSF/Data/LP_PlasmaCSF_230226.xlsx")
names(bf)[names(bf) == "PL_pTau217_pgml_Lumipulse"] <- "ptau217"
bf$ptau217 <- as.numeric(as.character(bf$ptau217))
names(bf)[names(bf) == "csf_clinical_routine_Abeta42_40_ratio_x10"] <- "csf_4240"
bf$csf_4240 <- as.numeric(as.character(bf$csf_4240))
bf$csf_4240 <- bf$csf_4240/10
bf$PL_217_42 <- (bf$ptau217)/(bf$PL_Ab42_pgml_Lumipulse)
bf <- bf %>% drop_na(csf_4240, ptau217)
bf$petstat_CL <- as.factor(ifelse(bf$CL_fnc_ber_com_composite >= 24, 1, 0)) 

###BF-MC
df.val <- bf %>% filter(bf$Study.y %in% "VALIDATE")

###BF-PC
df.ad <- bf %>% filter(bf$Study.y %in% "ADetect")

##ADNI
adni <- read_xlsx("~/Documents/Projects/LP vs CSF/Data/ADNI_VR2.xlsx")
adni$PL_217_42 <- (adni$ptau217)/(adni$PL_AB42) #PL_AB42_c
adni$CSF_18142 <- ifelse(adni$ABETA42 > 1700, (adni$PTAU_csf)/1700, (adni$PTAU_csf)/adni$ABETA42)
adni$cog_stat <- ifelse(adni$CDRSB_bl >= 0.5, 1, 0)
adni <- adni %>% filter(cog_stat == 1) 
adni$petstat_CL <- as.factor(ifelse(adni$CENTILOIDS >= 24, 1, 0)) 




----------------------- #1. Comparing regulatory approved cutoffs in each cohort ----------------------

#1A. AUCs----

auc_calc  <- function(df, pet_status, csf_marker, output_xlsx, output_pdf, marker_colors, marker_labels) {
  
  biomarkers <- c("ptau217","PL_217_42",csf_marker)
  df$petstat <- as.factor(df[[pet_status]])
  
  #DeLong comparisons
  results_delong <- list()
  
  for (i in seq_along(biomarkers)){
    for (j in seq_along(biomarkers)){ 
      if (i < j) {
        predictor1 <- biomarkers[i]
        predictor2 <- biomarkers[j]
        
        formula1 <- as.formula(paste("petstat~",predictor1, collapse=""))
        formula2 <- as.formula(paste("petstat~",predictor2, collapse=""))
        
        roc_curve1 <- pROC::roc(formula1, data=df, quiet=T)
        roc_curve2 <- pROC::roc(formula2, data=df,quiet=T)
        
        res <- roc.test(roc_curve1, roc_curve2, method = "delong")
        
        long <- data.frame(
          comp = paste(predictor1, predictor2, sep=" vs. "),
          auc1a = round(roc_curve1$auc,3),
          auc2a = round(roc_curve2$auc,3),
          pval = round(res$p.value,3),
          pval_cor = p.adjust(res$p.value, method = "fdr"),
          auc1_ci_low = round(ci.auc(roc_curve1)[1],3),
          auc1_ci_high =round(ci.auc(roc_curve1)[3],3),
          auc2_ci_low =round(ci.auc(roc_curve2)[1],3),
          auc2_ci_high = round(ci.auc(roc_curve2)[3],3)
        )
        
        results_delong[[length(results_delong) + 1]] <- long  # <- this was missing
        
      }
    }
  }
  results_delong <- do.call(rbind.data.frame, results_delong)
  results_delong$auc_order <- c(results_delong$auc1a[1],results_delong$auc1a[3], results_delong$auc2a[3])
  results_delong$auc_order_ci_low <- c(results_delong$auc1_ci_low[1],results_delong$auc1_ci_low[3], results_delong$auc2_ci_low[3])
  results_delong$auc_order_ci_high <- c(results_delong$auc1_ci_high[1],results_delong$auc1_ci_high[3], results_delong$auc2_ci_high[3])
  
  openxlsx::write.xlsx(results_delong,
                       file=output_xlsx)
  
  #AUC Plot
  roc_curves <- lapply(biomarkers, function(b) {
    pROC::roc(df$petstat ~ df[[b]], quiet = T)
  })
  
  p <- ggroc(roc_curves) +
    scale_color_manual(values = marker_colors,
                       labels = marker_labels, guide="none") +
    labs(title = element_blank(), x = "Specificity", y = "Sensitivity", color = "Assay") +
    scale_x_reverse()+
    scale_y_continuous()+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
    theme_classic()+
    theme(panel.border = element_blank(),
          axis.text.y = element_text(size = 5,color="grey30"),
          axis.text.x = element_text(size = 5,color="grey30"),
          axis.title.x = element_text(size=6,color="black"),
          axis.title.y = element_text(size=6),color="black")
  
  y_positions <- seq(0.25, by = -0.1, length.out = length(roc_curves))
  for(k in seq_along(roc_curves)){
    ci_k <- ci.auc(roc_curves[[k]])
    p <- p + annotate(
      "text", x = .3, y = y_positions[k], 
      label = sprintf("%.2f (%.2f, %.2f)", auc(roc_curves[[k]]), ci_k[1], ci_k[3]),
      size = 2, color = marker_colors[k]
    )
  }
  
  ggsave(output_pdf, plot = p, device = "pdf",width = 35, dpi=500, height = 35, units = "mm")
}

#Run per cohort

auc_calc(
  df = df.val, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/DeLong_BFMC_VR.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/AUC_BFMC_VR.pdf",
  marker_colors = c("#4ca481", "#c52851", "#295dbf"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF AB42/40")
)

auc_calc(
  df = df.ad, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/DeLong_BFPC_VR.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/AUC_BFPC_VR.pdf",
  marker_colors = c("#4ca481", "#c52851", "#295dbf"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF AB42/40")
)

auc_calc(
  df = adni, 
  pet_status = "VR",
  csf_marker = "CSF_18142", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/DeLong_ADNI_VR.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/AUC_ADNI_VR.pdf",
  marker_colors = c("#4ca481", "#c52851", "#4ebbd6"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF p-tau181/AB42")
)



#1B. Diagnostic accuracy data, comparisons----

dta_calc <- function(df, pet_status, csf_marker, csf_pos_cutoff, csf_pos2_cutoff = NULL, csf_direction = "lower", output_xlsx, output_pval_xlsx, output_table_xlsx){
  
  two_cutpoint <- !is.null(csf_pos2_cutoff)  #don't run if no 2nd cutoff available
  
  ##Set-up 
  
  df$petstat <- df[[pet_status]]
  
  df$ptau217ab42_pos <- ifelse(df$PL_217_42 >= 0.00738, 1, 0) 
  df$ptau217ab42_neg <- ifelse(df$PL_217_42 <= 0.00370, 0, 1)
  
  df$ptau217_pos <- ifelse(df$ptau217 >= 0.278, 1, 0) 
  df$ptau217_neg <- ifelse(df$ptau217 <= 0.161, 0, 1)
  df$ptau217_pos1 <- ifelse(df$ptau217 <= 0.201, 0, 1)
  
  
  if(csf_direction == "lower"){
    df$csf_pos <- ifelse(df[[csf_marker]] <= csf_pos_cutoff, 1, 0)
  }
  else {
    df$csf_pos <- ifelse(df[[csf_marker]] > csf_pos_cutoff, 1, 0)
  }
  
  
  if(two_cutpoint){
    if(csf_direction == "lower"){
      df$csf_pos2 <- ifelse(df[[csf_marker]] <= csf_pos2_cutoff, 1, 0)
    }
    else {
      df$csf_pos2 <- ifelse(df[[csf_marker]] > csf_pos2_cutoff, 1, 0)
    }
  }
  
  ##Function to calculate DTA
  compute_metrics <- function(true, pred) {
    true <- as.numeric(as.character(true))  # factor -> numeric
    pred <- as.numeric(as.character(pred))
    
    v_tp = sum(true == 1 & pred == 1)
    v_fp = sum(true == 0 & pred == 1)
    v_tn = sum(true == 0 & pred == 0)
    v_fn = sum(true == 1 & pred == 0)
    c(
      sens = v_tp / (v_tp + v_fn),
      spec = v_tn / (v_tn + v_fp),
      acc = (v_tp + v_tn) / (v_tp + v_tn + v_fn + v_fp),
      ppv = cutpointr::ppv(tp = v_tp, fp = v_fp, tn = v_tn, fn = v_fn),
      npv = cutpointr::npv(tp = v_tp, fp = v_fp, tn = v_tn, fn = v_fn),
      tp = v_tp, fp = v_fp, tn = v_tn, fn = v_fn)
  }
  
  
  #             SINGLE CUTOFF
  
  ##Calculate Accuracy, PPV, NPV and 95% CIs
  
  f_comp_stats <- function(data, indices) {
    d <- as.data.frame(data[indices,])
    
    r1 <- compute_metrics(d$petstat, d$ptau217_pos1)
    r2 <- compute_metrics(d$petstat, d$ptau217ab42_pos)
    r3 <- compute_metrics(d$petstat, d$csf_pos)
    
    #metrics & differences between them
    c(r1["acc"], r2["acc"], r3["acc"], #3
      r1["ppv"], r2["ppv"], r3["ppv"], #6
      r1["npv"], r2["npv"], r3["npv"], #9
      r1["sens"], r2["sens"], r3["sens"], #12
      r1["spec"], r2["spec"], r3["spec"], #15
      r1["tp"], r2["tp"], r3["tp"], #18
      r1["fp"], r2["fp"], r3["fp"], #21
      r1["tn"], r2["tn"], r3["tn"], #24
      r1["fn"], r2["fn"], r3["fn"], #27
      
      #difs
      r1["acc"]-r2["acc"], r1["acc"]-r3["acc"], r2["acc"]-r3["acc"],
      r1["ppv"]-r2["ppv"], r1["ppv"]-r3["ppv"], r2["ppv"]-r3["ppv"],      
      r1["npv"]-r2["npv"], r1["npv"]-r3["npv"], r2["npv"]-r3["npv"])
  }   
  
  set.seed(12345)
  boot_results <- boot(data = df, statistic = f_comp_stats, R = 2000)
  
  ##Calculate pvalues 1 cutoff
  dif_indices_1cp <- list(
    "Accuracy 217 vs. 217_AB42" = 28,
    "PPV pval 217 vs. 217_AB42" = 31,
    "NPV pval 217 vs. 217_AB42" = 34,
    "Accuracy 217 vs. CSF" = 29,
    "PPV pval 217 vs. CSF" = 32,
    "NPV pval 217 vs. CSF" = 35,
    "Accuracy 217_AB42 vs. CSF" = 30,
    "PPV pval 217_AB42 vs. CSF" = 33,
    "NPV pval 217_AB42 vs. CSF" = 36)
  
  pvals_1cp <- data.frame(
    Pval = numeric(length(dif_indices_1cp)),
    row.names = c("Accuracy 217 vs. 217_AB42","Accuracy 217 vs. CSF","Accuracy 217_AB42 vs. CSF",
                  "PPV pval 217 vs. 217_AB42","PPV pval 217_AB42 vs. CSF","PPV pval 217 vs. CSF",
                  "NPV pval 217 vs. 217_AB42", "NPV pval 217 vs. CSF",
                  "NPV pval 217_AB42 vs. CSF"))
  
  bootstrap_replicates <- boot_results$t
  for (stat_name in names(dif_indices_1cp)) {
    index <- dif_indices_1cp[[stat_name]]
    boot_diff <- boot_results$t0[index]
    results_underH0_diff <- bootstrap_replicates[, index] - mean(bootstrap_replicates[, index])
    p_value <- mean(abs(results_underH0_diff) >= abs(boot_diff))
    pvals_1cp[paste0(stat_name), "Pval"] <- round(p_value, 3)
    pvals_1cp$Pval_cor <- p.adjust(pvals_1cp$Pval, method = "fdr")
  }
  
  ##Add outcomes to dataframe (function to establish dataframe as desired)
  create_df_1cp <- function(method_name, boot_acc_index, boot_ppv_index, boot_npv_index, boot_sen_index, boot_spec_index,boot_tp_index,boot_fp_index,boot_tn_index,boot_fn_index) {
    data.frame(
      method = method_name,
      accuracy = round(boot_results$t0[boot_acc_index], 3),
      acc_ci_lower = round(quantile(boot_results$t[, boot_acc_index], c(0.025, 0.975), na.rm=T)[1], 3),
      acc_ci_upper = round(quantile(boot_results$t[, boot_acc_index], c(0.025, 0.975), na.rm=T)[2], 3),
      PPV = round(boot_results$t0[boot_ppv_index], 3),
      ppv_ci_lower = round(quantile(boot_results$t[, boot_ppv_index], c(0.025, 0.975), na.rm=T)[1], 3),
      ppv_ci_upper = round(quantile(boot_results$t[, boot_ppv_index], c(0.025, 0.975), na.rm=T)[2], 3),
      NPV = round(boot_results$t0[boot_npv_index], 3),
      npv_ci_lower = round(quantile(boot_results$t[, boot_npv_index], c(0.025, 0.975), na.rm=T)[1], 3),
      npv_ci_upper = round(quantile(boot_results$t[, boot_npv_index], c(0.025, 0.975), na.rm=T)[2], 3),
      Sensitivity =round(boot_results$t0[boot_sen_index], 3),
      sen_ci_lower = round(quantile(boot_results$t[, boot_sen_index], c(0.025, 0.975), na.rm=T)[1], 3),
      sen_ci_upper = round(quantile(boot_results$t[, boot_sen_index], c(0.025, 0.975), na.rm=T)[2], 3),
      Specificity = round(boot_results$t0[boot_spec_index], 3),
      spec_ci_lower = round(quantile(boot_results$t[, boot_spec_index], c(0.025, 0.975), na.rm=T)[1], 3),
      spec_ci_upper = round(quantile(boot_results$t[, boot_spec_index], c(0.025, 0.975), na.rm=T)[2], 3),
      TP = round(boot_results$t0[boot_tp_index], 0),
      FP = round(boot_results$t0[boot_fp_index], 0),
      TN = round(boot_results$t0[boot_tn_index], 0),
      FN = round(boot_results$t0[boot_fn_index], 0),
      intermediate_n_percentage = NA,
      intermediate_n_percentage_low = NA,
      intermediate_n_percentage_high = NA,
      stringsAsFactors = FALSE,
      row.names = "1 Cutpoint"
    )
  }
  
  
  results_1cp_1 <- create_df_1cp(method_name = "Plasma p-tau217",1,4,7,10,13,16,19,22,25)
  results_1cp_2 <- create_df_1cp(method_name = "Plasma p-tau217/AB42",2,5,8,11,14,17,20,23,26)
  results_1cp_3 <- create_df_1cp(method_name = "CSF", 3,6,9,12,15,18,21,24,27)
  merged_results_1cp <- rbind(results_1cp_1,results_1cp_2,results_1cp_3)
  
  
  #             TWO CUTOFFS
  
  pvals_intermediate <- NULL
  merged_2cp         <- NULL
  
  
  if(two_cutpoint){
    
    df <- df %>% mutate(
      greyzone_ptau217 = case_when(
        ptau217_neg == 0 ~ 1, ptau217_pos == 1 ~ 3, TRUE ~ 2),
      greyzone_ptau217ab42 = case_when(
        ptau217ab42_neg == 0 ~ 1, ptau217ab42_pos == 1 ~ 3, TRUE ~ 2),
      greyzone_csf = case_when(
        csf_pos == 0 ~ 1, csf_pos2 == 1 ~ 3, TRUE ~ 2)
    )
    
    
    ##Function to bootstrap %intermediate
    f_greyzone <- function(data, indices){
      d <- data[indices, ]
      
      d <- d %>% mutate(
        greyzone_ptau217 = case_when(
          ptau217_neg == 0 ~ 1, ptau217_pos == 1 ~ 3, TRUE ~ 2),
        greyzone_ptau217ab42 = case_when(
          ptau217ab42_neg == 0 ~ 1, ptau217ab42_pos == 1 ~ 3, TRUE ~ 2),
        greyzone_csf = case_when(
          csf_pos == 0 ~ 1, csf_pos2 == 1 ~ 3, TRUE ~ 2)
      )
      
      int_1 <- length(d$greyzone_ptau217 == 2)
      int_perc1 <- round((length(which(d$greyzone_ptau217 == 2)) / nrow(d)) * 100,3)
      
      int_2 <- length(d$greyzone_ptau217ab42 == 2)
      int_perc2 <- round((length(which(d$greyzone_ptau217ab42 == 2)) / nrow(d)) * 100,3)
      
      int_3 <- length(d$greyzone_csf == 2)
      int_perc3 <- round((length(which(d$greyzone_csf == 2)) / nrow(d)) * 100,3)
      
      return(c(int_perc1, int_perc2, int_perc3, 
               int_perc1 - int_perc2,int_perc1 - int_perc3, int_perc2 - int_perc3))
    }
    
    set.seed(12345)
    boot_grey <- boot(data = df, statistic = f_greyzone, R = 2000)
    
    ##Calculate pvalue %intermediate differences
    stat_indices_grey <- list(
      "Int% pval p-tau217/AB42 vs. p-tau217" = 4 , 
      "Int% pval p-tau217/AB42 vs. CSF" = 5, 
      "Int% pval p-tau217 vs.CSF" = 6)
    
    pvals_intermediate <- data.frame(
      Pval = numeric(length(stat_indices_grey)),
      row.names = c("Int% pval p-tau217/AB42 vs. p-tau217","Int% pval p-tau217/AB42 vs. CSF","Int% pval p-tau217 vs.CSF"))
    
    bootstrap_replicates_grey <- boot_grey$t
    
    for (stat_name in names(stat_indices_grey)) {
      index <- stat_indices_grey[[stat_name]]
      boot_diff <- boot_grey$t0[index]
      results_underH0_diff <- bootstrap_replicates_grey[, index] - mean(bootstrap_replicates_grey[, index])
      p_value <- mean(abs(results_underH0_diff) >= abs(boot_diff))
      pvals_intermediate[paste0(stat_name), "Pval"] <- round(p_value, 3)
      pvals_intermediate$Pval_cor <- p.adjust(pvals_intermediate$Pval, method = "fdr")
    }
    
    ##Calculate Accuracy, PPV, NPV and 95% CIs
    
    f_comp_stats_2cp <- function(data, indices) {
      d <- as.data.frame(data[indices,])
      
      d <- d %>% mutate(
        greyzone_ptau217 = case_when(
          ptau217_neg == 0 ~ 1, ptau217_pos == 1 ~ 3, TRUE ~ 2),
        greyzone_ptau217ab42 = case_when(
          ptau217ab42_neg == 0 ~ 1, ptau217ab42_pos == 1 ~ 3, TRUE ~ 2),
        greyzone_csf = case_when(
          csf_pos == 0 ~ 1, csf_pos2 == 1 ~ 3, TRUE ~ 2)
      )
      
      compute_2cp <- function(d, grey_col, pred_col){
        dx <- d %>% filter(!!sym(grey_col) != 2)
        compute_metrics(dx$petstat, dx[[pred_col]])
      }
      
      r1 <- compute_2cp(d, "greyzone_ptau217", "ptau217_pos")
      r2 <- compute_2cp(d, "greyzone_ptau217ab42", "ptau217ab42_pos")
      r3 <- compute_2cp(d, "greyzone_csf", "csf_pos")
      
      #metrics & differences between them
      c(r1["acc"], r2["acc"], r3["acc"], #3
        r1["ppv"], r2["ppv"], r3["ppv"], #6
        r1["npv"], r2["npv"], r3["npv"], #9
        r1["sens"], r2["sens"], r3["sens"], #12
        r1["spec"], r2["spec"], r3["spec"], #15
        r1["tp"], r2["tp"], r3["tp"], #18
        r1["fp"], r2["fp"], r3["fp"], #21
        r1["tn"], r2["tn"], r3["tn"], #24
        r1["fn"], r2["fn"], r3["fn"], #27
        
        #difs
        r1["acc"]-r2["acc"], r1["acc"]-r3["acc"], r2["acc"]-r3["acc"],
        r1["ppv"]-r2["ppv"], r1["ppv"]-r3["ppv"], r2["ppv"]-r3["ppv"],      
        r1["npv"]-r2["npv"], r1["npv"]-r3["npv"], r2["npv"]-r3["npv"])
    }   
    
    set.seed(12345)
    boot_results_2cp <- boot(data = df, statistic = f_comp_stats_2cp, R = 2000)
    
    ##Calculate pvalues 1 cutoff
    dif_indices_2cp <- list(
      "Accuracy 217 vs. 217_AB42" = 28,
      "PPV pval 217 vs. 217_AB42" = 31,
      "NPV pval 217 vs. 217_AB42" = 34,
      "Accuracy 217 vs. CSF" = 29,
      "PPV pval 217 vs. CSF" = 32,
      "NPV pval 217 vs. CSF" = 35,
      "Accuracy 217_AB42 vs. CSF" = 30,
      "PPV pval 217_AB42 vs. CSF" = 33,
      "NPV pval 217_AB42 vs. CSF" = 36)
    
    pvals_2cp <- data.frame(
      Pval = numeric(length(dif_indices_2cp)),
      row.names = c("Accuracy 217 vs. 217_AB42","Accuracy 217 vs. CSF","Accuracy 217_AB42 vs. CSF",
                    "PPV pval 217 vs. 217_AB42","PPV pval 217_AB42 vs. CSF","PPV pval 217 vs. CSF",
                    "NPV pval 217 vs. 217_AB42", "NPV pval 217 vs. CSF",
                    "NPV pval 217_AB42 vs. CSF"))
    
    bootstrap_replicates <- boot_results_2cp$t
    for (stat_name in names(dif_indices_2cp)) {
      index <- dif_indices_2cp[[stat_name]]
      boot_diff <- boot_results_2cp$t0[index]
      results_underH0_diff <- bootstrap_replicates[, index] - mean(bootstrap_replicates[, index])
      p_value <- mean(abs(results_underH0_diff) >= abs(boot_diff))
      pvals_2cp[paste0(stat_name), "Pval"] <- round(p_value, 3)
      pvals_2cp$Pval_cor <- p.adjust(pvals_2cp$Pval, method = "fdr")
    }
    
    ##Add outcomes to dataframe (function to establish dataframe as desired)
    create_df_2cp <- function(method_name, boot_acc_index, boot_ppv_index, boot_npv_index, boot_sen_index, boot_spec_index,boot_tp_index,boot_fp_index,boot_tn_index,boot_fn_index, grey_index) {
      data.frame(
        method = method_name,
        accuracy = round(boot_results_2cp$t0[boot_acc_index], 3),
        acc_ci_lower = round(quantile(boot_results_2cp$t[, boot_acc_index], c(0.025, 0.975), na.rm=T)[1], 3),
        acc_ci_upper = round(quantile(boot_results_2cp$t[, boot_acc_index], c(0.025, 0.975), na.rm=T)[2], 3),
        PPV = round(boot_results_2cp$t0[boot_ppv_index], 3),
        ppv_ci_lower = round(quantile(boot_results_2cp$t[, boot_ppv_index], c(0.025, 0.975), na.rm=T)[1], 3),
        ppv_ci_upper = round(quantile(boot_results_2cp$t[, boot_ppv_index], c(0.025, 0.975), na.rm=T)[2], 3),
        NPV = round(boot_results_2cp$t0[boot_npv_index], 3),
        npv_ci_lower = round(quantile(boot_results_2cp$t[, boot_npv_index], c(0.025, 0.975), na.rm=T)[1], 3),
        npv_ci_upper = round(quantile(boot_results_2cp$t[, boot_npv_index], c(0.025, 0.975), na.rm=T)[2], 3),
        Sensitivity =round(boot_results_2cp$t0[boot_sen_index], 3),
        sen_ci_lower = round(quantile(boot_results_2cp$t[, boot_sen_index], c(0.025, 0.975), na.rm=T)[1], 3),
        sen_ci_upper = round(quantile(boot_results_2cp$t[, boot_sen_index], c(0.025, 0.975), na.rm=T)[2], 3),
        Specificity = round(boot_results_2cp$t0[boot_spec_index], 3),
        spec_ci_lower = round(quantile(boot_results_2cp$t[, boot_spec_index], c(0.025, 0.975), na.rm=T)[1], 3),
        spec_ci_upper = round(quantile(boot_results_2cp$t[, boot_spec_index], c(0.025, 0.975), na.rm=T)[2], 3),
        TP = round(boot_results_2cp$t0[boot_tp_index], 0),
        FP = round(boot_results_2cp$t0[boot_fp_index], 0),
        TN = round(boot_results_2cp$t0[boot_tn_index], 0),
        FN = round(boot_results_2cp$t0[boot_fn_index], 0),
        intermediate_n_percentage = round(boot_grey$t0[grey_index], 3),
        intermediate_n_percentage_low = round(quantile(boot_grey$t[, grey_index], c(0.025, 0.975), na.rm=T)[1], 3),
        intermediate_n_percentage_high = round(quantile(boot_grey$t[, grey_index], c(0.025, 0.975), na.rm=T)[2], 3),
        stringsAsFactors = FALSE,
        row.names = "1 Cutpoint"
      )
    }
    
    results_1cp_1 <- create_df_2cp(method_name = "Plasma p-tau217",1,4,7,10,13,16,19,22,25,1)
    results_1cp_2 <- create_df_2cp(method_name = "Plasma p-tau217/AB42",2,5,8,11,14,17,20,23,26,2)
    results_1cp_3 <- create_df_2cp(method_name = "CSF", 3,6,9,12,15,18,21,24,27,3)
    merged_results_2cp <- rbind(results_1cp_1, results_1cp_2, results_1cp_3)
  }
  
  #             COMBINE
  
  
  if(two_cutpoint){
    merged_comps <- rbind(merged_results_1cp, merged_results_2cp)
    merged_comps$Cutpoints <- c("1 Cutpoint","1 Cutpoint","1 Cutpoint",
                                "2 Cutpoints","2 Cutpoints","2 Cutpoints")
    pvalues_all <- rbind(pvals_1cp, pvals_intermediate, pvals_2cp)
  } else {
    merged_comps <- merged_results_1cp
    merged_comps$Cutpoints <- "1 Cutpoint"
    pvalues_all <- pvals_1cp
    
  }
  
  openxlsx::write.xlsx(merged_comps,file=output_xlsx, rowNames =T)
  openxlsx::write.xlsx(pvalues_all,file=output_pval_xlsx, rowNames =T)
  
  table_df <- data.frame(
    Measure = merged_comps$method,
    Accuracy = sprintf("%.2f (%.2f-%.2f)", merged_comps$accuracy, merged_comps$acc_ci_lower, merged_comps$acc_ci_upper),
    PPV = sprintf("%.2f (%.2f-%.2f)", merged_comps$PPV, merged_comps$ppv_ci_lower, merged_comps$ppv_ci_upper),
    NPV = sprintf("%.2f (%.2f-%.2f)", merged_comps$NPV, merged_comps$npv_ci_lower, merged_comps$npv_ci_upper),
    Sensitivity = sprintf("%.2f (%.2f-%.2f)", merged_comps$Sensitivity, merged_comps$sen_ci_lower, merged_comps$sen_ci_upper),
    Specificity =sprintf("%.2f (%.2f-%.2f)", merged_comps$Specificity, merged_comps$spec_ci_lower, merged_comps$spec_ci_upper), 
    Gre= sprintf("%.2f (%.2f-%.2f)", merged_comps$intermediate_n_percentage, merged_comps$intermediate_n_percentage_low, merged_comps$intermediate_n_percentage_high),
    TP =sprintf("%.2f", merged_comps$TP),
    FP =sprintf("%.2f", merged_comps$FP),
    TN =sprintf("%.2f", merged_comps$TN), 
    FN =sprintf("%.2f", merged_comps$FN))
  openxlsx::write.xlsx(table_df,file=output_table_xlsx, rowNames =T)
}


#Run per dataset 

dta_calc(
  df = df.val, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR_table.xlsx"
)

dta_calc(
  df = df.ad, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR_table.xlsx"
)

dta_calc(
  df = adni, 
  pet_status = "VR",
  csf_marker = "CSF_18142", 
  csf_pos_cutoff = 0.028,
  csf_pos2_cutoff = NULL,
  csf_direction = "upper",
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_VR.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_VR_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_VR_table.xlsx"
)


#1C. Main figures----

##1CP plot----
create_plot_1cp <- function(df, csf_color, output_pdf){
  
  make_plot_data <- function(rows, cols, measure_name) {
    d <- cbind(df[rows, "method"], df[rows, cols])
    colnames(d) <- c("Method", "Measure1", "ci_low", "ci_high")
    d$Measure <- measure_name
    d
  }
  
  plotdata <- rbind(
    make_plot_data(1:3, c("accuracy", "acc_ci_lower", "acc_ci_upper"), "Accuracy"),
    make_plot_data(1:3, c("PPV", "ppv_ci_lower", "ppv_ci_upper"), "PPV"),
    make_plot_data(1:3, c("NPV", "npv_ci_lower", "npv_ci_upper"), "NPV")
  )
  
  plotdata$Method2 <- as.factor(c("1", "2", "3"))
  
  
  group_order <- c("NPV", "PPV", "Accuracy")
  method_order <- c("3", "2", "1")
  tick_positions <- seq(40, 100, by = 10)
  colors <- c("1" = "#4ca481", "2" = "#c52851", "3" = csf_color)
  
  
  #Forest plot Acc, PPV, NPV
  p1 <- ggplot(data=plotdata, aes(y = factor(Measure, levels=group_order), x = Measure1 * 100, 
                                  color=factor(Method2, levels=method_order), fill=factor(Method2, levels=method_order)))+
    geom_vline(xintercept = tick_positions, color = "grey80", linewidth = 0.5, alpha=0.4) +
    geom_errorbar(aes(xmin = ci_low * 100, xmax = ci_high * 100),
                  width = 0, linewidth = 1.05, position = position_dodge(width = 0.6)) +
    geom_point(position = position_dodge(width = 0.6), size = 2, stroke = 0.7)+
    scale_x_continuous(limits=c(50,110),breaks = seq(50,100, by=10),expand =c(0,0))+
    labs(x = "Percentage", y = "") +
    geom_text(aes(label = sprintf("%.f (%.f-%.f)", Measure1 * 100, ci_low * 100, ci_high * 100),
                  group=factor(Method2, levels=method_order), x=106), position = position_dodge(width = 0.6), size = 2) +
    scale_color_manual(values=colors, guide="none")+
    scale_fill_manual(values=colors, guide="none")+
    theme_classic() + 
    theme(axis.text.x = element_text(size = 7,color="grey30"),
          axis.text.y = element_text(size = 7, color="black"),
          axis.title.x = element_text(size = 7,color="black"),
          axis.ticks.y = element_blank(),
          legend.text = element_text(size = 6,color="black"),  
          legend.title = element_blank(), 
          legend.position = "bottom")
  
  ggsave(output_pdf, p1, device = "pdf",width = 60, dpi=500, height = 52, units = "mm")
}



create_plot_1cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR.xlsx"),
                csf_color = "#295dbf", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate//Test/Res_1cp.pdf")

create_plot_1cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR.xlsx"),
                csf_color = "#295dbf", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_1cp.pdf")

create_plot_1cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_VR.xlsx"),
                csf_color = "#4ebbd6", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_1cp.pdf")


##2CP plot----
create_plot_2cp <- function(df, csf_color, output_pdf, output_pdf_grey){
  
  make_plot_data <- function(rows, cols, measure_name) {
    d <- cbind(df[rows, "method"], df[rows, cols])
    colnames(d) <- c("Method", "Measure1", "ci_low", "ci_high")
    d$Measure <- measure_name
    d
  }
  
  plotdata <- rbind(
    make_plot_data(4:6, c("accuracy", "acc_ci_lower", "acc_ci_upper"), "Accuracy"),
    make_plot_data(4:6, c("PPV", "ppv_ci_lower", "ppv_ci_upper"), "PPV"),
    make_plot_data(4:6, c("NPV", "npv_ci_lower", "npv_ci_upper"), "NPV")
  )
  
  plotdata$Method2 <- as.factor(c("1", "2", "3"))
  
  
  group_order <- c("NPV", "PPV", "Accuracy")
  method_order <- c("3", "2", "1")
  tick_positions <- seq(40, 100, by = 10)
  colors <- c("1" = "#4ca481", "2" = "#c52851", "3" = csf_color)
  
  
  #Forest plot Acc, PPV, NPV
  p1 <- ggplot(data=plotdata, aes(y = factor(Measure, levels=group_order), x = Measure1 * 100, 
                                  color=factor(Method2, levels=method_order), fill=factor(Method2, levels=method_order)))+
    geom_vline(xintercept = tick_positions, color = "grey80", linewidth = 0.5, alpha=0.4) +
    geom_errorbar(aes(xmin = ci_low * 100, xmax = ci_high * 100),
                  width = 0, linewidth = 1.05, position = position_dodge(width = 0.6)) +
    geom_point(position = position_dodge(width = 0.6), size = 2, stroke = 0.7)+
    scale_x_continuous(limits=c(50,110),breaks = seq(50,100, by=10),expand =c(0,0))+
    labs(x = "Percentage", y = "") +
    geom_text(aes(label = sprintf("%.f (%.f-%.f)", Measure1 * 100, ci_low * 100, ci_high * 100),
                  group=factor(Method2, levels=method_order), x=106), position = position_dodge(width = 0.6), size = 2) +
    scale_color_manual(values=colors, guide="none")+
    scale_fill_manual(values=colors, guide="none")+
    theme_classic() + 
    theme(axis.text.x = element_text(size = 7,color="grey30"),
          axis.text.y = element_text(size = 7, color="black"),
          axis.title.x = element_text(size = 7,color="black"),
          axis.ticks.y = element_blank(),
          legend.text = element_text(size = 6,color="black"),  
          legend.title = element_blank(), 
          legend.position = "bottom")
  
  ggsave(output_pdf, p1, device = "pdf",width = 60, dpi=500, height = 52, units = "mm")
  
  
  #Intermediate plot
  plot_grey <- data.frame(
    Method = df$method[4:6],
    Measure1 = df$intermediate_n_percentage[4:6],
    ci_low = df$intermediate_n_percentage_low[4:6],
    ci_high = df$intermediate_n_percentage_high[4:6]
  )
  
  method_levels <- plot_grey$Method
  colors <- c("Plasma p-tau217" = "#4ca481", "Plasma p-tau217/AB42" = "#c52851", "CSF" = csf_color)
  
  
  grey <- ggplot(data=plot_grey, aes(x=factor(Method,levels = method_levels),y=Measure1, 
                                     fill=factor(Method, levels = method_levels)))+
    geom_bar(color="black",stat = "identity", alpha=0.8, width=0.6,  position = position_dodge(width = 0.1),)+
    geom_errorbar(aes(x=factor(Method,level = method_levels), 
                      ymax = ci_high, ymin = ci_low),width=0.2, size =0.5,position=position_dodge(0.1),  color="black")+
    scale_y_continuous(limits=c(0,50),breaks = seq(0,50, by=10),expand =c(0,0))+
    geom_text(aes(label = sprintf("%.1f", Measure1)), position = position_dodge(width=0.75),
              vjust=-3, size = 2) +
    scale_fill_manual(values = colors, guide="none") +
    labs(title = element_blank(),x = "",y = "Intermediate values (%)") +
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=6),
          axis.title.x = element_text(size=7),
          axis.title.y= element_text(size=7),
          panel.border = element_rect(linewidth = 0.5, fill = NA, color=NA))  
  ggsave(output_pdf_grey, grey, device = "pdf",width = 40, dpi=500, height = 40, units = "mm")
}


create_plot_2cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR.xlsx"),
                csf_color = "#295dbf", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate//Test/Res_1cp.pdf",
                output_pdf_grey = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_grey.pdf")

create_plot_2cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR.xlsx"),
                csf_color = "#295dbf", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_1cp.pdf",
                output_pdf_grey = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_grey.pdf")


----------------------- #2. Comparing within-cohort cutoffs in each cohort ----------------------

#2B. Diagnostic accuracy data, comparisons within----

dta_calc_within <- function(df, pet_status, 
                            ptau217_spec90_cutoff, ptau217_spec95_cutoff, ptau217_sens95_cutoff,
                            ptau217ab42_spec90_cutoff, ptau217ab42_spec95_cutoff, ptau217ab42_sens95_cutoff,
                            csf_spec90_cutoff, csf_spec95_cutoff, csf_sens95_cutoff,
                            csf_marker, csf_pos_cutoff, csf_direction = "lower", 
                            output_xlsx, output_pval_xlsx, output_table_xlsx){
  
  
  ##Set-up 
  two_cutpoint <- !is.null(csf_spec95_cutoff)  #don't run if no 2nd cutoff available
  
  df$petstat <- df[[pet_status]]
  
  # p-tau217
  df$ptau217_sens95 <- ifelse(df$ptau217 > ptau217_sens95_cutoff, 1, 0)
  df$ptau217_spec90 <- ifelse(df$ptau217 > ptau217_spec90_cutoff, 1, 0)
  df$ptau217_spec95 <- ifelse(df$ptau217 > ptau217_spec95_cutoff, 1, 0)
  
  # p-tau217/AB42
  df$ptau217ab42_sens95 <- ifelse(df$PL_217_42 > ptau217ab42_sens95_cutoff, 1, 0)
  df$ptau217ab42_spec90 <- ifelse(df$PL_217_42 > ptau217ab42_spec90_cutoff, 1, 0)
  df$ptau217ab42_spec95 <- ifelse(df$PL_217_42 > ptau217ab42_spec95_cutoff, 1, 0)
  
  
  df$csf_spec90 <- ifelse(df$csf_4240 < 0.055, 1, 0) 
  df$csf_spec95 <- ifelse(df$csf_4240 < 0.051, 1, 0)
  df$csf_sens95 <- ifelse(df$csf_4240 < 0.066, 1, 0)
  
  # CSF direction varies per cohort
  if(csf_direction == "lower") {
    df$csf_sens95 <- ifelse(df[[csf_marker]] < csf_sens95_cutoff, 1, 0)
    df$csf_spec90 <- ifelse(df[[csf_marker]] < csf_spec90_cutoff, 1, 0)
    df$csf_spec95 <- ifelse(df[[csf_marker]] < csf_spec95_cutoff, 1, 0)
  } else {
    df$csf_sens95 <- ifelse(df[[csf_marker]] > csf_sens95_cutoff, 1, 0)
    df$csf_spec90 <- ifelse(df[[csf_marker]] > csf_spec90_cutoff, 1, 0)
    df$csf_spec95 <- ifelse(df[[csf_marker]] > csf_spec95_cutoff, 1, 0)
  }
  
  
  
  ##Function to calculate DTA
  compute_metrics <- function(true, pred) {
    true <- as.numeric(as.character(true))  # factor -> numeric
    pred <- as.numeric(as.character(pred))
    
    v_tp = sum(true == 1 & pred == 1)
    v_fp = sum(true == 0 & pred == 1)
    v_tn = sum(true == 0 & pred == 0)
    v_fn = sum(true == 1 & pred == 0)
    c(
      sens = v_tp / (v_tp + v_fn),
      spec = v_tn / (v_tn + v_fp),
      acc = (v_tp + v_tn) / (v_tp + v_tn + v_fn + v_fp),
      ppv = cutpointr::ppv(tp = v_tp, fp = v_fp, tn = v_tn, fn = v_fn),
      npv = cutpointr::npv(tp = v_tp, fp = v_fp, tn = v_tn, fn = v_fn),
      tp = v_tp, fp = v_fp, tn = v_tn, fn = v_fn)
  }
  
  
  ###### SINGLE CUTOFF #########
  
  
  ##Calculate Accuracy, PPV, NPV and 95% CIs
  
  f_comp_stats <- function(data, indices) {
    d <- as.data.frame(data[indices,])
    
    r1 <- compute_metrics(d$petstat, d$ptau217_spec90)
    r2 <- compute_metrics(d$petstat, d$ptau217ab42_spec90)
    r3 <- compute_metrics(d$petstat, d$csf_spec90)
    
    #metrics & differences between them
    c(r1["acc"], r2["acc"], r3["acc"], #3
      r1["ppv"], r2["ppv"], r3["ppv"], #6
      r1["npv"], r2["npv"], r3["npv"], #9
      r1["sens"], r2["sens"], r3["sens"], #12
      r1["spec"], r2["spec"], r3["spec"], #15
      r1["tp"], r2["tp"], r3["tp"], #18
      r1["fp"], r2["fp"], r3["fp"], #21
      r1["tn"], r2["tn"], r3["tn"], #24
      r1["fn"], r2["fn"], r3["fn"], #27
      
      #difs
      r1["acc"]-r2["acc"], r1["acc"]-r3["acc"], r2["acc"]-r3["acc"],
      r1["ppv"]-r2["ppv"], r1["ppv"]-r3["ppv"], r2["ppv"]-r3["ppv"],      
      r1["npv"]-r2["npv"], r1["npv"]-r3["npv"], r2["npv"]-r3["npv"])
  }   
  
  set.seed(12345)
  boot_results <- boot(data = df, statistic = f_comp_stats, R = 2000)
  
  ##Calculate pvalues 1 cutoff
  dif_indices_1cp <- list(
    "Accuracy 217 vs. 217_AB42" = 28,
    "PPV pval 217 vs. 217_AB42" = 31,
    "NPV pval 217 vs. 217_AB42" = 34,
    "Accuracy 217 vs. CSF" = 29,
    "PPV pval 217 vs. CSF" = 32,
    "NPV pval 217 vs. CSF" = 35,
    "Accuracy 217_AB42 vs. CSF" = 30,
    "PPV pval 217_AB42 vs. CSF" = 33,
    "NPV pval 217_AB42 vs. CSF" = 36)
  
  pvals_1cp <- data.frame(
    Pval = numeric(length(dif_indices_1cp)),
    row.names = c("Accuracy 217 vs. 217_AB42","Accuracy 217 vs. CSF","Accuracy 217_AB42 vs. CSF",
                  "PPV pval 217 vs. 217_AB42","PPV pval 217_AB42 vs. CSF","PPV pval 217 vs. CSF",
                  "NPV pval 217 vs. 217_AB42", "NPV pval 217 vs. CSF",
                  "NPV pval 217_AB42 vs. CSF"))
  
  bootstrap_replicates <- boot_results$t
  for (stat_name in names(dif_indices_1cp)) {
    index <- dif_indices_1cp[[stat_name]]
    boot_diff <- boot_results$t0[index]
    results_underH0_diff <- bootstrap_replicates[, index] - mean(bootstrap_replicates[, index])
    p_value <- mean(abs(results_underH0_diff) >= abs(boot_diff))
    pvals_1cp[paste0(stat_name), "Pval"] <- round(p_value, 3)
    pvals_1cp$Pval_cor <- p.adjust(pvals_1cp$Pval, method = "fdr")
  }
  
  ##Add outcomes to dataframe (function to establish dataframe as desired)
  create_df_1cp <- function(method_name, boot_acc_index, boot_ppv_index, boot_npv_index, boot_sen_index, boot_spec_index,boot_tp_index,boot_fp_index,boot_tn_index,boot_fn_index) {
    data.frame(
      method = method_name,
      accuracy = round(boot_results$t0[boot_acc_index], 3),
      acc_ci_lower = round(quantile(boot_results$t[, boot_acc_index], c(0.025, 0.975), na.rm=T)[1], 3),
      acc_ci_upper = round(quantile(boot_results$t[, boot_acc_index], c(0.025, 0.975), na.rm=T)[2], 3),
      PPV = round(boot_results$t0[boot_ppv_index], 3),
      ppv_ci_lower = round(quantile(boot_results$t[, boot_ppv_index], c(0.025, 0.975), na.rm=T)[1], 3),
      ppv_ci_upper = round(quantile(boot_results$t[, boot_ppv_index], c(0.025, 0.975), na.rm=T)[2], 3),
      NPV = round(boot_results$t0[boot_npv_index], 3),
      npv_ci_lower = round(quantile(boot_results$t[, boot_npv_index], c(0.025, 0.975), na.rm=T)[1], 3),
      npv_ci_upper = round(quantile(boot_results$t[, boot_npv_index], c(0.025, 0.975), na.rm=T)[2], 3),
      Sensitivity =round(boot_results$t0[boot_sen_index], 3),
      sen_ci_lower = round(quantile(boot_results$t[, boot_sen_index], c(0.025, 0.975), na.rm=T)[1], 3),
      sen_ci_upper = round(quantile(boot_results$t[, boot_sen_index], c(0.025, 0.975), na.rm=T)[2], 3),
      Specificity = round(boot_results$t0[boot_spec_index], 3),
      spec_ci_lower = round(quantile(boot_results$t[, boot_spec_index], c(0.025, 0.975), na.rm=T)[1], 3),
      spec_ci_upper = round(quantile(boot_results$t[, boot_spec_index], c(0.025, 0.975), na.rm=T)[2], 3),
      TP = round(boot_results$t0[boot_tp_index], 0),
      FP = round(boot_results$t0[boot_fp_index], 0),
      TN = round(boot_results$t0[boot_tn_index], 0),
      FN = round(boot_results$t0[boot_fn_index], 0),
      intermediate_n_percentage = NA,
      intermediate_n_percentage_low = NA,
      intermediate_n_percentage_high = NA,
      stringsAsFactors = FALSE,
      row.names = "1 Cutpoint"
    )
  }
  
  
  results_1cp_1 <- create_df_1cp(method_name = "Plasma p-tau217",1,4,7,10,13,16,19,22,25)
  results_1cp_2 <- create_df_1cp(method_name = "Plasma p-tau217/AB42",2,5,8,11,14,17,20,23,26)
  results_1cp_3 <- create_df_1cp(method_name = "CSF", 3,6,9,12,15,18,21,24,27)
  merged_results_1cp <- rbind(results_1cp_1,results_1cp_2,results_1cp_3)
  
  
  #####TWO CUTOFFS#####
  
  pvals_intermediate <- NULL
  merged_2cp         <- NULL
  
  
  if(two_cutpoint){
    
    df <- df %>% mutate(
      greyzone_ptau217 = case_when(
        ptau217_sens95 == 0 ~ 1, ptau217_spec95 == 1 ~ 3, TRUE ~ 2),
      greyzone_ptau217ab42 = case_when(
        ptau217ab42_sens95 == 0 ~ 1, ptau217ab42_spec95 == 1 ~ 3, TRUE ~ 2),
      greyzone_csf = case_when(
        csf_sens95 == 0 ~ 1, csf_spec95 == 1 ~ 3, TRUE ~ 2)
    )
    
    
    ##Function to bootstrap %intermediate
    f_greyzone <- function(data, indices){
      d <- data[indices, ]
      
      d <- d %>% mutate(
        greyzone_ptau217 = case_when(
          ptau217_sens95 == 0 ~ 1, ptau217_spec95 == 1 ~ 3, TRUE ~ 2),
        greyzone_ptau217ab42 = case_when(
          ptau217ab42_sens95 == 0 ~ 1, ptau217ab42_spec95 == 1 ~ 3, TRUE ~ 2),
        greyzone_csf = case_when(
          csf_sens95 == 0 ~ 1, csf_spec95 == 1 ~ 3, TRUE ~ 2)
      )
      
      int_1 <- length(d$greyzone_ptau217 == 2)
      int_perc1 <- round((length(which(d$greyzone_ptau217 == 2)) / nrow(d)) * 100,3)
      
      int_2 <- length(d$greyzone_ptau217ab42 == 2)
      int_perc2 <- round((length(which(d$greyzone_ptau217ab42 == 2)) / nrow(d)) * 100,3)
      
      int_3 <- length(d$greyzone_csf == 2)
      int_perc3 <- round((length(which(d$greyzone_csf == 2)) / nrow(d)) * 100,3)
      
      return(c(int_perc1, int_perc2, int_perc3, 
               int_perc1 - int_perc2,int_perc1 - int_perc3, int_perc2 - int_perc3))
    }
    
    set.seed(12345)
    boot_grey <- boot(data = df, statistic = f_greyzone, R = 2000)
    
    ##Calculate pvalue %intermediate differences
    stat_indices_grey <- list(
      "Int% pval p-tau217/AB42 vs. p-tau217" = 4 , 
      "Int% pval p-tau217/AB42 vs. CSF" = 5, 
      "Int% pval p-tau217 vs.CSF" = 6)
    
    pvals_intermediate <- data.frame(
      Pval = numeric(length(stat_indices_grey)),
      row.names = c("Int% pval p-tau217/AB42 vs. p-tau217","Int% pval p-tau217/AB42 vs. CSF","Int% pval p-tau217 vs.CSF"))
    
    bootstrap_replicates_grey <- boot_grey$t
    
    for (stat_name in names(stat_indices_grey)) {
      index <- stat_indices_grey[[stat_name]]
      boot_diff <- boot_grey$t0[index]
      results_underH0_diff <- bootstrap_replicates_grey[, index] - mean(bootstrap_replicates_grey[, index])
      p_value <- mean(abs(results_underH0_diff) >= abs(boot_diff))
      pvals_intermediate[paste0(stat_name), "Pval"] <- round(p_value, 3)
      pvals_intermediate$Pval_cor <- p.adjust(pvals_intermediate$Pval, method = "fdr")
    }
    
    ##Calculate Accuracy, PPV, NPV and 95% CIs
    
    f_comp_stats_2cp <- function(data, indices) {
      d <- as.data.frame(data[indices,])
      
      d <- d %>% mutate(
        greyzone_ptau217 = case_when(
          ptau217_sens95 == 0 ~ 1, ptau217_spec95 == 1 ~ 3, TRUE ~ 2),
        greyzone_ptau217ab42 = case_when(
          ptau217ab42_sens95 == 0 ~ 1, ptau217ab42_spec95 == 1 ~ 3, TRUE ~ 2),
        greyzone_csf = case_when(
          csf_sens95 == 0 ~ 1, csf_spec95 == 1 ~ 3, TRUE ~ 2)
      )
      
      compute_2cp <- function(d, grey_col, pred_col){
        dx <- d %>% filter(!!sym(grey_col) != 2)
        compute_metrics(dx$petstat, dx[[pred_col]])
      }
      
      r1 <- compute_2cp(d, "greyzone_ptau217", "ptau217_spec90")
      r2 <- compute_2cp(d, "greyzone_ptau217ab42", "ptau217ab42_spec90")
      r3 <- compute_2cp(d, "greyzone_csf", "csf_spec90")
      
      #metrics & differences between them
      c(r1["acc"], r2["acc"], r3["acc"], #3
        r1["ppv"], r2["ppv"], r3["ppv"], #6
        r1["npv"], r2["npv"], r3["npv"], #9
        r1["sens"], r2["sens"], r3["sens"], #12
        r1["spec"], r2["spec"], r3["spec"], #15
        r1["tp"], r2["tp"], r3["tp"], #18
        r1["fp"], r2["fp"], r3["fp"], #21
        r1["tn"], r2["tn"], r3["tn"], #24
        r1["fn"], r2["fn"], r3["fn"], #27
        
        #difs
        r1["acc"]-r2["acc"], r1["acc"]-r3["acc"], r2["acc"]-r3["acc"],
        r1["ppv"]-r2["ppv"], r1["ppv"]-r3["ppv"], r2["ppv"]-r3["ppv"],      
        r1["npv"]-r2["npv"], r1["npv"]-r3["npv"], r2["npv"]-r3["npv"])
    }   
    
    set.seed(12345)
    boot_results_2cp <- boot(data = df, statistic = f_comp_stats_2cp, R = 2000)
    
    ##Calculate pvalues 1 cutoff
    dif_indices_2cp <- list(
      "Accuracy 217 vs. 217_AB42" = 28,
      "PPV pval 217 vs. 217_AB42" = 31,
      "NPV pval 217 vs. 217_AB42" = 34,
      "Accuracy 217 vs. CSF" = 29,
      "PPV pval 217 vs. CSF" = 32,
      "NPV pval 217 vs. CSF" = 35,
      "Accuracy 217_AB42 vs. CSF" = 30,
      "PPV pval 217_AB42 vs. CSF" = 33,
      "NPV pval 217_AB42 vs. CSF" = 36)
    
    pvals_2cp <- data.frame(
      Pval = numeric(length(dif_indices_2cp)),
      row.names = c("Accuracy 217 vs. 217_AB42","Accuracy 217 vs. CSF","Accuracy 217_AB42 vs. CSF",
                    "PPV pval 217 vs. 217_AB42","PPV pval 217_AB42 vs. CSF","PPV pval 217 vs. CSF",
                    "NPV pval 217 vs. 217_AB42", "NPV pval 217 vs. CSF",
                    "NPV pval 217_AB42 vs. CSF"))
    
    bootstrap_replicates <- boot_results_2cp$t
    for (stat_name in names(dif_indices_2cp)) {
      index <- dif_indices_2cp[[stat_name]]
      boot_diff <- boot_results_2cp$t0[index]
      results_underH0_diff <- bootstrap_replicates[, index] - mean(bootstrap_replicates[, index])
      p_value <- mean(abs(results_underH0_diff) >= abs(boot_diff))
      pvals_2cp[paste0(stat_name), "Pval"] <- round(p_value, 3)
      pvals_2cp$Pval_cor <- p.adjust(pvals_2cp$Pval, method = "fdr")
    }
    
    ##Add outcomes to dataframe (function to establish dataframe as desired)
    create_df_2cp <- function(method_name, boot_acc_index, boot_ppv_index, boot_npv_index, boot_sen_index, boot_spec_index,boot_tp_index,boot_fp_index,boot_tn_index,boot_fn_index, grey_index) {
      data.frame(
        method = method_name,
        accuracy = round(boot_results_2cp$t0[boot_acc_index], 3),
        acc_ci_lower = round(quantile(boot_results_2cp$t[, boot_acc_index], c(0.025, 0.975), na.rm=T)[1], 3),
        acc_ci_upper = round(quantile(boot_results_2cp$t[, boot_acc_index], c(0.025, 0.975), na.rm=T)[2], 3),
        PPV = round(boot_results_2cp$t0[boot_ppv_index], 3),
        ppv_ci_lower = round(quantile(boot_results_2cp$t[, boot_ppv_index], c(0.025, 0.975), na.rm=T)[1], 3),
        ppv_ci_upper = round(quantile(boot_results_2cp$t[, boot_ppv_index], c(0.025, 0.975), na.rm=T)[2], 3),
        NPV = round(boot_results_2cp$t0[boot_npv_index], 3),
        npv_ci_lower = round(quantile(boot_results_2cp$t[, boot_npv_index], c(0.025, 0.975), na.rm=T)[1], 3),
        npv_ci_upper = round(quantile(boot_results_2cp$t[, boot_npv_index], c(0.025, 0.975), na.rm=T)[2], 3),
        Sensitivity =round(boot_results_2cp$t0[boot_sen_index], 3),
        sen_ci_lower = round(quantile(boot_results_2cp$t[, boot_sen_index], c(0.025, 0.975), na.rm=T)[1], 3),
        sen_ci_upper = round(quantile(boot_results_2cp$t[, boot_sen_index], c(0.025, 0.975), na.rm=T)[2], 3),
        Specificity = round(boot_results_2cp$t0[boot_spec_index], 3),
        spec_ci_lower = round(quantile(boot_results_2cp$t[, boot_spec_index], c(0.025, 0.975), na.rm=T)[1], 3),
        spec_ci_upper = round(quantile(boot_results_2cp$t[, boot_spec_index], c(0.025, 0.975), na.rm=T)[2], 3),
        TP = round(boot_results_2cp$t0[boot_tp_index], 0),
        FP = round(boot_results_2cp$t0[boot_fp_index], 0),
        TN = round(boot_results_2cp$t0[boot_tn_index], 0),
        FN = round(boot_results_2cp$t0[boot_fn_index], 0),
        intermediate_n_percentage = round(boot_grey$t0[grey_index], 3),
        intermediate_n_percentage_low = round(quantile(boot_grey$t[, grey_index], c(0.025, 0.975), na.rm=T)[1], 3),
        intermediate_n_percentage_high = round(quantile(boot_grey$t[, grey_index], c(0.025, 0.975), na.rm=T)[2], 3),
        stringsAsFactors = FALSE,
        row.names = "1 Cutpoint"
      )
    }
    
    results_1cp_1 <- create_df_2cp(method_name = "Plasma p-tau217",1,4,7,10,13,16,19,22,25,1)
    results_1cp_2 <- create_df_2cp(method_name = "Plasma p-tau217/AB42",2,5,8,11,14,17,20,23,26,2)
    results_1cp_3 <- create_df_2cp(method_name = "CSF", 3,6,9,12,15,18,21,24,27,3)
    merged_results_2cp <- rbind(results_1cp_1, results_1cp_2, results_1cp_3)
  }
  
  
  ###### COMBINE ########
  
  
  if(two_cutpoint){
    merged_comps <- rbind(merged_results_1cp, merged_results_2cp)
    merged_comps$Cutpoints <- c("1 Cutpoint","1 Cutpoint","1 Cutpoint",
                                "2 Cutpoints","2 Cutpoints","2 Cutpoints")
    pvalues_all <- rbind(pvals_1cp, pvals_intermediate, pvals_2cp)
  } else {
    merged_comps <- merged_results_1cp
    merged_comps$Cutpoints <- "1 Cutpoint"
    pvalues_all <- pvals_1cp
    
  }
  
  openxlsx::write.xlsx(merged_comps,file=output_xlsx, rowNames =T)
  openxlsx::write.xlsx(pvalues_all,file=output_pval_xlsx, rowNames =T)
  
  table_df <- data.frame(
    Measure = merged_comps$method,
    Accuracy = sprintf("%.2f (%.2f-%.2f)", merged_comps$accuracy, merged_comps$acc_ci_lower, merged_comps$acc_ci_upper),
    PPV = sprintf("%.2f (%.2f-%.2f)", merged_comps$PPV, merged_comps$ppv_ci_lower, merged_comps$ppv_ci_upper),
    NPV = sprintf("%.2f (%.2f-%.2f)", merged_comps$NPV, merged_comps$npv_ci_lower, merged_comps$npv_ci_upper),
    Sensitivity = sprintf("%.2f (%.2f-%.2f)", merged_comps$Sensitivity, merged_comps$sen_ci_lower, merged_comps$sen_ci_upper),
    Specificity =sprintf("%.2f (%.2f-%.2f)", merged_comps$Specificity, merged_comps$spec_ci_lower, merged_comps$spec_ci_upper), 
    Gre= sprintf("%.2f (%.2f-%.2f)", merged_comps$intermediate_n_percentage, merged_comps$intermediate_n_percentage_low, merged_comps$intermediate_n_percentage_high),
    TP =sprintf("%.2f", merged_comps$TP),
    FP =sprintf("%.2f", merged_comps$FP),
    TN =sprintf("%.2f", merged_comps$TN), 
    FN =sprintf("%.2f", merged_comps$FN))
  openxlsx::write.xlsx(table_df,file=output_table_xlsx, rowNames =T)
}


#Run per dataset 
dta_calc_within(
  df = df.val,
  pet_status = "VR_overall",
  ptau217_spec90_cutoff = 0.253, ptau217_spec95_cutoff = 0.310, ptau217_sens95_cutoff  = 0.190, 
  ptau217ab42_spec90_cutoff = 0.009, ptau217ab42_spec95_cutoff = 0.011,ptau217ab42_sens95_cutoff = 0.007,
  csf_spec90_cutoff = 0.055, csf_spec95_cutoff = 0.051, csf_sens95_cutoff = 0.066, csf_marker = "csf_4240", csf_direction = "lower",
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR_Within.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR_pval_Within.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR_table_Within.xlsx")

dta_calc_within(
  df = df.ad,
  pet_status = "VR_overall",
  ptau217_spec90_cutoff = 0.253, ptau217_spec95_cutoff = 0.310, ptau217_sens95_cutoff  = 0.190, 
  ptau217ab42_spec90_cutoff = 0.009, ptau217ab42_spec95_cutoff = 0.011,ptau217ab42_sens95_cutoff = 0.007,
  csf_spec90_cutoff = 0.055, csf_spec95_cutoff = 0.051, csf_sens95_cutoff = 0.066, csf_marker = "csf_4240", csf_direction = "lower",
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR_Within.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR_pval_Within.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR_table_Within.xlsx")

dta_calc_within(
  df = adni,
  pet_status = "VR",
  ptau217_spec90_cutoff = 0.17, ptau217_spec95_cutoff = 0.20, ptau217_sens95_cutoff  = 0.09, 
  ptau217ab42_spec90_cutoff = 0.007, ptau217ab42_spec95_cutoff = 0.009,ptau217ab42_sens95_cutoff = 0.003,
  csf_spec90_cutoff = 0.019, csf_spec95_cutoff = 0.022, csf_sens95_cutoff = 0.021, csf_marker = "CSF_18142", csf_direction = "upper",
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_VR_Within.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_VR_pval_Within.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_VR_table_Within.xlsx")


----------------------- #3. CENTILOID status ----------------------

#3A. AUC----
auc_calc(
  df = df.val, 
  pet_status = "petstat_CL",
  csf_marker = "csf_4240", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/DeLong_BFMC_CL.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/AUC_BFMC_CL.pdf",
  marker_colors = c("#4ca481", "#c52851", "#295dbf"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF AB42/40")
)

auc_calc(
  df = df.ad, 
  pet_status = "petstat_CL",
  csf_marker = "csf_4240", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/DeLong_BFPC_CL.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/AUC_BFPC_CL.pdf",
  marker_colors = c("#4ca481", "#c52851", "#295dbf"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF AB42/40")
)

auc_calc(
  df = adni, 
  pet_status = "petstat_CL",
  csf_marker = "CSF_18142", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/DeLong_ADNI_CL.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/AUC_ADNI_CL.pdf",
  marker_colors = c("#4ca481", "#c52851", "#4ebbd6"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF p-tau181/AB42")
)

#3B. Diagnostic accuracy CL ----

dta_calc(
  df = df.val, 
  pet_status = "petstat_CL",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_CL.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_CL_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_CL_table.xlsx"
)

dta_calc(
  df = df.ad, 
  pet_status = "petstat_CL",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_CL.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_CL_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_CL_table.xlsx"
)

dta_calc(
  df = adni, 
  pet_status = "petstat_CL",
  csf_marker = "CSF_18142", 
  csf_pos_cutoff = 0.028,
  csf_pos2_cutoff = NULL,
  csf_direction = "upper",
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_CL.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_CL_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_CL_table.xlsx"
)

#3C. Plots CL----

##1cp
create_plot_1cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_CL.xlsx"),
                csf_color = "#295dbf", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate//Test/Res_1cp_CL.pdf")

create_plot_1cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_CL.xlsx"),
                csf_color = "#295dbf", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_1cp_CL.pdf")

create_plot_1cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_CL.xlsx"),
                csf_color = "#4ebbd6", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Test/Res_1cp_CL.pdf")

##2cp

create_plot_2cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_CL.xlsx"),
                csf_color = "#295dbf", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate//Test/Res_1cp_CL.pdf",
                output_pdf_grey = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_grey_CL.pdf")

create_plot_2cp(df = read_xlsx("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_CL.xlsx"),
                csf_color = "#295dbf", 
                output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_1cp_CL.pdf",
                output_pdf_grey = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_grey_CL.pdf")




----------------------- #4. CI only ----------------------
#4A. AUC----

df.val_ci <- df.val %>% filter(!cognitive_status_baseline_variable %in% c("Normal", "SCD"))

auc_calc(
  df = df.val_ci, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/DeLong_BFMC_VR_CI.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/AUC_BFMC_VR_CI.pdf",
  marker_colors = c("#4ca481", "#c52851", "#295dbf"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF AB42/40")
)

df.ad_ci <- df.ad %>% filter(!cognitive_status_baseline_variable %in% c("Normal", "SCD"))

auc_calc(
  df = df.ad_ci, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/DeLong_BFPC_VR_CI.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/AUC_BFPC_VR_CI.pdf",
  marker_colors = c("#4ca481", "#c52851", "#295dbf"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF AB42/40")
)


#4B. Diagnostic accuracy CI ----

dta_calc(
  df = df.val_ci, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR_CI.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR_CI_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Test/Res_VR_CI_table.xlsx"
)

dta_calc(
  df = df.ad_ci, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR_CI.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR_CI_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Test/Res_VR_CI_table.xlsx"
)


#5. SCD only ----------------------

df.val_scd <- df.val %>% filter(cognitive_status_baseline_variable %in% c("SCD"))
df.ad_scd <- df.ad %>% filter(cognitive_status_baseline_variable %in% c("SCD"))

df.pool_scd <- rbind(df.ad_scd, df.val_scd)

#5A. AUC----
auc_calc(
  df = df.pool_scd, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/DeLong_BFMC_VR_SCD.xlsx",
  output_pdf = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/AUC_BFMC_VR_SCD.pdf",
  marker_colors = c("#4ca481", "#c52851", "#295dbf"),
  marker_labels = c("Plasma p-tau217", "Plasma p-tau217/AB42", "CSF AB42/40")
)


#5B. Diagnostic accuracy SCD ----

dta_calc(
  df = df.pool_scd, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_SCD.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_SCD_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_SCD_table.xlsx"
)


#6. Age strat ----------------------

df.pool_bf <- rbind(df.ad, df.ad)
df.pool_bf$age_category <- ifelse(df.pool_bf$age_blood_test < median(df.pool_bf$age_blood_test), 0, 1)
table(df.pool_bf$age_category)
median(na.omit(df$age_blood_test))


###AGE 1
df_age1 <- df.pool_bf %>% filter(age_category == 0)

dta_calc(
  df = df_age1, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_age1.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_age1_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_age1_table.xlsx"
)

###AGE 2
df_age2 <- df.pool_bf %>% filter(age_category == 1)

dta_calc(
  df = df_age2, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_age2.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_age2_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_age2_table.xlsx"
)

#7. Sex strat ----------------------

df.pool_bf <- rbind(df.ad, df.ad)

###Male
df_male <- df.pool_bf %>% filter(gender_baseline_variable == 0)

dta_calc(
  df = df_male, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_M.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_M_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_M_table.xlsx"
)

###Female
df_female <- df.pool_bf %>% filter(gender_baseline_variable == 1)

dta_calc(
  df = df_female, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_F.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_F_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_F_table.xlsx"
)




#8. Outliers ----------------------

df.pool_bf <- rbind(df.ad, df.ad)

###Male
df_male <- df.pool_bf %>% filter(gender_baseline_variable == 0)

dta_calc(
  df = df_male, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_M.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_M_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_M_table.xlsx"
)

###Female
df_female <- df.pool_bf %>% filter(gender_baseline_variable == 1)

dta_calc(
  df = df_female, 
  pet_status = "VR_overall",
  csf_marker = "csf_4240", 
  csf_pos_cutoff = 0.072,
  csf_pos2_cutoff = 0.058,
  output_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_F.xlsx",
  output_pval_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_F_pval.xlsx",
  output_table_xlsx = "~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Test/Res_VR_F_table.xlsx"
)






