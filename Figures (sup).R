#SUPPLEMENTARY FIGURES

#Load libraries

library(readxl)
library(tidyverse)
library(pROC)
library(boot)
library(gtsummary)
library(ggpubr)
library(ggbeeswarm)
library(pROC)
library(see)

#BioFINDER-Memory Clinic----
df <- read_xlsx("~/Documents/Projects/LP vs CSF/Data/LP_PlasmaCSF_300625.xlsx")
df.val <- df %>% filter(df$Study %in% "VALIDATE")
df <- df.val #124

#change names etc.
names(df)[names(df) == "Plasma_Ptau217P_pgmL_Lumi"] <- "ptau217"
df$ptau217 <- as.numeric(as.character(df$ptau217))
names(df)[names(df) == "csf_clinical_routine_Abeta42_40_ratio_x10"] <- "csf_4240"
df$csf_4240 <- as.numeric(as.character(df$csf_4240))
df$csf_4240 <- df$csf_4240/10
#pass passing-bablok transformation for Lumipulse AB42 assay
df$Plasma_AB142_pgmL_Lumi_c <- 1.119 + 1.114*df$Plasma_AB142_pgmL_Lumi
df$PL_217_42 <- (df$ptau217)/(df$Plasma_AB142_pgmL_Lumi_c)

df <- df %>% select(sid_bf2, sid_diag, VR_overall, CL_fnc_ber_com_composite, csf_4240, ptau217, PL_217_42, Plasma_AB142_pgmL_Lumi_c, Plasma_AB142_pgmL_Lumi)
df <- df %>% drop_na(csf_4240, ptau217)

##S1 concordance biomarkers----

df$petstat <- as.factor(df$VR_overall)

cor.test(df$ptau217, df$PL_217_42, method="spearman")
##p-tau217 and p-tau217/AB42
ggplot(df, aes(x = PL_217_42, y = ptau217, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_vline(xintercept=0.00738, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_hline(yintercept=0.23211, linetype="dotted", color="grey20", linewidth = 0.4)+
  ylab("Plasma p-tau217 (pg/mL)")+
  xlab("Plasma p-tau217/AB42")+
  scale_x_continuous(breaks = round(seq(0, 0.1, by = 0.05), 3), limits = c(0, 0.1)) +
  scale_y_continuous(breaks = seq(0, 2.5, by = .5), limits = c(0, 2.5),  labels = function(x) sprintf("%.1f", x)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"), guide = "none") + #003f5c
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 0.075, y = 1.75, label = "rho = 0.97, p = < 2.2e-16", size = 2)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Conc_ptau217_ab42.pdf",device = "pdf",width = 50, dpi=500, height = 45, units = "mm")

##p-tau217 and CSF

cor.test(df$ptau217, df$csf_4240, method="spearman")

ggplot(df, aes(x = ptau217, y = csf_4240, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.072, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_vline(xintercept=0.23211, linetype="dotted", color="grey20", linewidth = 0.4)+
  xlab("Plasma p-tau217 (pg/mL)")+
  ylab("CSF AB42/40")+
  scale_y_continuous(breaks = round(seq(0.025, 0.125, by = 0.05), 3), limits = c(0.025, 0.125)) +
  scale_x_continuous(breaks = seq(0, 2, by = .5), limits = c(0, 2),  labels = function(x) sprintf("%.1f", x)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"),guide="none")+
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 1.75, y = 0.100, label = "rho = -0.64, p = < 2.2e-16", size = 2)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Conc_ptau217_CSF.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")

##p-tau217/AB42 and CSF

cor.test(df$PL_217_42, df$csf_4240, method="spearman")

ggplot(df, aes(x = PL_217_42, y = csf_4240, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.072, linetype="dotted", color="grey20", linewidth = 0.6)+
  geom_vline(xintercept=0.00738, linetype="dotted", color="grey20", linewidth = 0.6)+
  xlab("Plasma p-tau217/AB42")+
  ylab("CSF AB42/40")+
  scale_y_continuous(breaks = round(seq(0.025, 0.125, by = 0.05), 3), limits = c(0.025, 0.125)) +
  scale_x_continuous(breaks = round(seq(0, 0.08, by = 0.02), 3), limits = c(0, 0.02)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"),guide="none")+
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 0.015, y = 0.08, label = "rho = -0.66, p = < 2.2e-16", size = 2)
ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Conc_ptau217ab42_CSF.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")




##S8 correlations biomarkers and centiloids----

ggplot(df, aes(x = CL_fnc_ber_com_composite, y = ptau217)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#4ca481") +
  geom_smooth(method = "lm", color = "#4ca481", fill = "#4ca481") +
  scale_x_continuous(breaks = seq(-25, 175, by = 50), limits = c(-25, 175)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1), limits = c(0, 3.5)) +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  labs(x = "AB-PET (Centiloids)", y = "Plasma p-tau217") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=3, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/LPptau217_Cent_correlations.pdf",device = "pdf",width = 50,
       dpi=500, height = 40, units = "mm")

ggplot(df, aes(x = CL_fnc_ber_com_composite, y = PL_217_42)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#c52851") +
  geom_smooth(method = "lm", color = "#c52851", fill = "#c52851") +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  scale_x_continuous(breaks = seq(-25, 175, by = 50), limits = c(-25, 175)) +
  scale_y_continuous(breaks = seq(0, .1, by = .05), limits = c(0, .122)) +
  labs(x = "AB-PET (Centiloids)", y = "Plasma  p-tau217/AB42") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=.075, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/LPptau217AB42_Cent_correlations.pdf",device = "pdf",
       width = 50, dpi=500, height = 40, units = "mm")

ggplot(df, aes(x = CL_fnc_ber_com_composite, y = csf_4240)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#295dbf") +
  scale_x_continuous(breaks = seq(-25, 175, by = 50), limits = c(-25, 175)) +
  scale_y_continuous(breaks = seq(0.025, .125, by = .025), limits = c(0.025, .125)) +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  labs(x = "AB-PET (Centiloids)", y = "CSF AB42/40") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=0.1, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Csf_Cent_correlations.pdf",device = "pdf",
       width = 50, dpi=500, height = 40, units = "mm")

##S9 correlations between VR and PET----
ggplot(df, aes(x = as.factor(VR_overall), y = CL_fnc_ber_com_composite, fill = as.factor(VR_overall)))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.08681, ymax = 0.23211,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  geom_hline(yintercept=24, linetype="dotted", color="grey20") +
  #geom_hline(yintercept=0.08681, linetype="dotted", color="grey20") +
  scale_y_continuous(breaks = round(seq(-25, 175, by = 25), 1), limits = c(-25, 185)) +
  scale_x_discrete(labels=c("0" ="VR-","1" ="VR+")) +
  labs(x = element_blank(), y = "AB-PET Centiloids") +
  scale_fill_manual(values = c("#afafaf", "#af639e"), guide = "none") +
  scale_color_manual(values = c("#afafaf", "#af639e"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 150)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Correlations_CL_VR.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")



#in BioFINDER-Primary Care----

df <- read_xlsx("~/Documents/Projects/LP vs CSF/Data/LP_PlasmaCSF_300625.xlsx")

#change names etc.
names(df)[names(df) == "Plasma_Ptau217P_pgmL_Lumi"] <- "ptau217"
df$ptau217 <- as.numeric(as.character(df$ptau217))
names(df)[names(df) == "csf_clinical_routine_Abeta42_40_ratio_x10"] <- "csf_4240"
df$csf_4240 <- as.numeric(as.character(df$csf_4240))
df$csf_4240 <- df$csf_4240/10
#PB transformation for AB42
df$Plasma_AB142_pgmL_Lumi_c <- 1.119 + 1.114*df$Plasma_AB142_pgmL_Lumi
df$PL_217_42 <- (df$ptau217)/(df$Plasma_AB142_pgmL_Lumi_c)

#select dataset
df.ad <- df %>% filter(df$Study %in% "ADetect")
df <- df.ad #124
df <- df %>% select(VR_overall, CL_fnc_ber_com_composite, csf_4240, ptau217, PL_217_42, Plasma_AB142_pgmL_Lumi_c, Plasma_AB142_pgmL_Lumi)
df <- df %>% drop_na(csf_4240,ptau217)

##S1 concordance biomarkers----

df$petstat <- as.factor(df$VR_overall)

cor.test(df$ptau217, df$PL_217_42, method="spearman")
##p-tau217 and p-tau217/AB42
ggplot(df, aes(x = PL_217_42, y = ptau217, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_vline(xintercept=0.00738, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_hline(yintercept=0.23211, linetype="dotted", color="grey20", linewidth = 0.4)+
  ylab("Plasma p-tau217 (pg/mL)")+
  xlab("Plasma p-tau217/AB42")+
  scale_x_continuous(breaks = round(seq(0, 0.15, by = 0.05), 3), limits = c(0, 0.15)) +
  scale_y_continuous(breaks = seq(0, 4, by = 1), limits = c(0, 4),  labels = function(x) sprintf("%.1f", x)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"), guide = "none") + #003f5c
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 0.075, y = 1.75, label = "rho = 0.98, p = < 2.2e-16", size = 2)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/FDA approach/Conc_ptau217_ab42.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")

##p-tau217 and CSF

cor.test(df$ptau217, df$csf_4240, method="spearman")

ggplot(df, aes(x = ptau217, y = csf_4240, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.072, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_vline(xintercept=0.23211, linetype="dotted", color="grey20", linewidth = 0.4)+
  xlab("Plasma p-tau217 (pg/mL)")+
  ylab("CSF AB42/40")+
  scale_y_continuous(breaks = round(seq(0.025, 0.125, by = 0.05), 3), limits = c(0.025, 0.125)) +
  scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(0, 4),  labels = function(x) sprintf("%.1f", x)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"),guide="none")+
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 1.75, y = 0.100, label = "rho = -0.64, p = < 2.2e-16", size = 2)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/FDA approach/Conc_ptau217_CSF.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")

##p-tau217/AB42 and CSF

cor.test(df$PL_217_42, df$csf_4240, method="spearman")

ggplot(df, aes(x = PL_217_42, y = csf_4240, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.072, linetype="dotted", color="grey20", linewidth = 0.6)+
  geom_vline(xintercept=0.00738, linetype="dotted", color="grey20", linewidth = 0.6)+
  xlab("Plasma p-tau217/AB42")+
  ylab("CSF AB42/40")+
  scale_y_continuous(breaks = round(seq(0.025, 0.125, by = 0.05), 3), limits = c(0.025, 0.125)) +
  scale_x_continuous(breaks = round(seq(0, 0.15, by = 0.05), 3), limits = c(0, 0.15)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"),guide="none")+
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 0.075, y = 0.08, label = "rho = -0.66, p = < 2.2e-16", size = 2)
ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/FDA approach/Conc_ptau217ab42_CSF.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")



##S8 correlations biomarkers and centiloids----

ggplot(df, aes(x = CL_fnc_ber_com_composite, y = ptau217)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#4ca481") +
  geom_smooth(method = "lm", color = "#4ca481", fill = "#4ca481") +
  scale_x_continuous(breaks = seq(-25, 175, by = 50), limits = c(-25, 185)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1), limits = c(0, 4)) +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  labs(x = "AB-PET (Centiloids)", y = "Plasma p-tau217") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=3, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/LPptau217_Cent_correlations.pdf",device = "pdf",
       width = 50, dpi=500, height = 40, units = "mm")

ggplot(df, aes(x = CL_fnc_ber_com_composite, y = PL_217_42)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#c52851") +
  geom_smooth(method = "lm", color = "#c52851", fill = "#c52851") +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  scale_x_continuous(breaks = seq(-25, 175, by = 50), limits = c(-25, 185)) +
  scale_y_continuous(breaks = seq(0, .1, by = .05), limits = c(0, .122)) +
  labs(x = "AB-PET (Centiloids)", y = "Plasma  p-tau217/AB42") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=.075, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/LPptau217AB42_Cent_correlations.pdf",device = "pdf",
       width = 50, dpi=500, height = 40, units = "mm")

ggplot(df, aes(x = CL_fnc_ber_com_composite, y = csf_4240)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#295dbf") +
  scale_x_continuous(breaks = seq(-25, 175, by = 50), limits = c(-25, 185)) +
  scale_y_continuous(breaks = seq(0.025, .125, by = .025), limits = c(0.025, .125)) +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  labs(x = "AB-PET (Centiloids)", y = "CSF AB42/40") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=0.1, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Csf_Cent_correlations.pdf",device = "pdf",
       width = 50, dpi=500, height = 40, units = "mm")

##S9 correlations between VR and PET----
df <- df %>% drop_na(CL_fnc_ber_com_composite)
ggplot(df, aes(x = as.factor(VR_overall), y = CL_fnc_ber_com_composite, fill = as.factor(VR_overall)))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.08681, ymax = 0.23211,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  geom_hline(yintercept=24, linetype="dotted", color="grey20") +
  #geom_hline(yintercept=0.08681, linetype="dotted", color="grey20") +
  scale_y_continuous(breaks = round(seq(-25, 175, by = 25), 1), limits = c(-25, 185)) +
  scale_x_discrete(labels=c("0" ="VR-","1" ="VR+")) +
  labs(x = element_blank(), y = "AB-PET Centiloids") +
  scale_fill_manual(values = c("#afafaf", "#af639e"), guide = "none") +
  scale_color_manual(values = c("#afafaf", "#af639e"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 150)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Primary care/Correlations_CL_VR.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")


##in ADNI----

df <- read_xlsx("~/Documents/Projects/LP vs CSF/Data/ADNI_VR2.xlsx")
df$PL_AB42_c <- 1.119 + 1.114*df$PL_AB42
#df$PL_217_42 <- (df$ptau217)/(df$Plasma_AB142_pgmL_Lumi)
df$PL_217_42 <- (df$ptau217)/(df$PL_AB42_c)
df$PL_217_42_n <- (df$ptau217)/(df$PL_AB42)
df$CSF_18142 <- ifelse(df$ABETA42 > 1700, (df$PTAU_csf)/1700, (df$PTAU_csf)/df$ABETA42)
df$cog_stat <- ifelse(df$CDRSB_bl >= 0.5, 1, 0)
df <- df %>% filter(cog_stat == 1) 


##S1 concordance biomarkers----

df$petstat <- as.factor(df$VR)

cor.test(df$ptau217, df$PL_217_42, method="spearman")

##p-tau217 and p-tau217/AB42
ggplot(df, aes(x = PL_217_42, y = ptau217, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_vline(xintercept=0.00738, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_hline(yintercept=0.23211, linetype="dotted", color="grey20", linewidth = 0.4)+
  ylab("Plasma p-tau217 (pg/mL)")+
  xlab("Plasma p-tau217/AB42")+
  scale_x_continuous(breaks = round(seq(0, 0.125, by = 0.05), 3), limits = c(0, 0.125)) +
  scale_y_continuous(breaks = seq(0, 2, by = .5), limits = c(0, 2),  labels = function(x) sprintf("%.1f", x)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"), guide = "none") + #003f5c
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 0.075, y = 1.75, label = "rho = 0.95, p = < 2.2e-16", size = 2)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/FDA approach/Conc_ptau217_ab42.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")

##p-tau217 and CSF

cor.test(df$ptau217, df$CSF_18142, method="spearman")

ggplot(df, aes(x = ptau217, y = CSF_18142, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.072, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_vline(xintercept=0.28, linetype="dotted", color="grey20", linewidth = 0.4)+
  xlab("Plasma p-tau217 (pg/mL)")+
  ylab("CSF p-tau181/AB42")+
  scale_y_continuous(breaks = round(seq(0.0, .2, by = 0.1), 3), limits = c(0.0, .2)) +
  scale_x_continuous(breaks = seq(0, 3.2, by = 1), limits = c(0, 3.2),  labels = function(x) sprintf("%.1f", x)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"),guide="none")+
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 1.75, y = 0.100, label = "rho = 0.75, p = < 2.2e-16", size = 2)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/FDA approach/Conc_ptau217_CSF.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")

##p-tau217/AB42 and CSF

cor.test(df$PL_217_42, df$CSF_18142, method="spearman")

ggplot(df, aes(x = PL_217_42, y = CSF_18142, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.028, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_vline(xintercept=0.00738, linetype="dotted", color="grey20", linewidth = 0.6)+
  xlab("Plasma p-tau217/AB42")+
  ylab("CSF p-tau181/AB42")+
  scale_y_continuous(breaks = round(seq(0.0, .2, by = 0.1), 3), limits = c(0.0, .2)) +
  scale_x_continuous(breaks = round(seq(0, .1, by = 0.05), 3), limits = c(0, .1)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"),guide="none")+
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 0.075, y = 0.08, label = "rho = 0.72, p = < 2.2e-16", size = 2)
ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/FDA approach/Conc_ptau217ab42_CSF.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")

##S8 correlations biomarkers and centiloids----

ggplot(df, aes(x = CENTILOIDS, y = ptau217)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#4ca481") +
  geom_smooth(method = "lm", color = "#4ca481", fill = "#4ca481") +
  scale_x_continuous(breaks = seq(-25, 175, by = 25), limits = c(-40, 175)) +
  scale_y_continuous(breaks = seq(0, 3, by = 1), limits = c(0, 3.1)) +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  labs(x = "AB-PET (Centiloids)", y = "Plasma p-tau217") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=3, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/LPptau217_Cent_correlations.pdf",device = "pdf",width = 50, dpi=500, height = 40, units = "mm")

ggplot(df, aes(x = CENTILOIDS, y = PL_217_42)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#c52851") +
  geom_smooth(method = "lm", color = "#c52851", fill = "#c52851") +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  scale_x_continuous(breaks = seq(-25, 175, by = 25), limits = c(-40, 175)) +
  scale_y_continuous(breaks = seq(0, .1, by = .05), limits = c(0, .122)) +
  labs(x = "AB-PET (Centiloids)", y = "Plasma  p-tau217/AB42") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=.075, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/LPptau217AB42_Cent_correlations.pdf",device = "pdf",width = 50, dpi=500, height = 40, units = "mm")

ggplot(df, aes(x = CENTILOIDS, y = CSF_18142)) +
  geom_point(size = 2, shape = 20, alpha = 0.85, color="#ea8a0f") +
  scale_x_continuous(breaks = seq(-25, 175, by = 50), limits = c(-25, 185)) +
  scale_y_continuous(breaks = seq(0.025, .2, by = .05), limits = c(0, .2)) +
  geom_vline(xintercept=24, linetype="dotted", color="black", alpha=0.9)+
  labs(x = "AB-PET (Centiloids)", y = "CSF p-tau181/AB42") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold", size=7), axis.title.y=element_text(face="bold", size=5), axis.text.x = element_text(size=6), axis.text.y=element_text(size=6))+
  stat_cor(method="spearman", label.x = 75, label.y=0.1, cor.coef.name = "rho", p.digits = 0.01)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Csf_Cent_correlations.pdf",device = "pdf",width = 50, dpi=500, height = 40, units = "mm")


##S9 correlations between VR and PET----
df <- df %>% drop_na(CENTILOIDS)

ggplot(df, aes(x = as.factor(VR), y = CENTILOIDS, fill = as.factor(VR)))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.08681, ymax = 0.23211,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  geom_hline(yintercept=24, linetype="dotted", color="grey20") +
  #geom_hline(yintercept=0.08681, linetype="dotted", color="grey20") +
  scale_y_continuous(breaks = round(seq(-25, 175, by = 25), 1), limits = c(-40, 175)) +
  scale_x_discrete(labels=c("0" ="VR-","1" ="VR+")) +
  labs(x = element_blank(), y = "AB-PET Centiloids") +
  scale_fill_manual(values = c("#afafaf", "#af639e"), guide = "none") +
  scale_color_manual(values = c("#afafaf", "#af639e"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 150)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/ADNI/Correlations_CL_VR.pdf",device = "pdf",
       width = 50, dpi=500, height = 45, units = "mm")
