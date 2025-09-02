#SUPPLEMENTARY FIGURES

#Load data----

##BioFINDER
df <- read_xlsx("~/Documents/Projects/LP vs CSF/Data/LP_PlasmaCSF_300625.xlsx")
#change names etc.
names(df)[names(df) == "Plasma_Ptau217P_pgmL_Lumi"] <- "ptau217"
df$ptau217 <- as.numeric(as.character(df$ptau217))
names(df)[names(df) == "csf_clinical_routine_Abeta42_40_ratio_x10"] <- "csf_4240"
df$csf_4240 <- as.numeric(as.character(df$csf_4240))
df$csf_4240 <- df$csf_4240/10
#pass passing-bablok transformation for Lumipulse AB42 assay
df$Plasma_AB142_pgmL_Lumi_c <- 1.119 + 1.114*df$Plasma_AB142_pgmL_Lumi
df$PL_217_42 <- (df$ptau217)/(df$Plasma_AB142_pgmL_Lumi_c)

###BF-Memory clinic
df.val <- df %>% filter(df$Study %in% "VALIDATE")
df.val <- df.val %>% select(sid_bf2, sid_diag, VR_overall, CL_fnc_ber_com_composite, csf_4240, ptau217, PL_217_42, Plasma_AB142_pgmL_Lumi_c, Plasma_AB142_pgmL_Lumi, petstat)
df.val <- df.val %>% drop_na(csf_4240, ptau217)

###BF-Primary care


##ADNI



df$petstat <- as.factor(df$VR_overall)


#S1 concordance biomarkers----

##in BFMC----

cor.test(df.val$ptau217, df.val$PL_217_42, method="spearman")
##p-tau217 and p-tau217/AB42
ggplot(df.val, aes(x = ptau217, y = PL_217_42, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.00738, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_vline(xintercept=0.23211, linetype="dotted", color="grey20", linewidth = 0.4)+
  xlab("Plasma p-tau217 (pg/mL)")+
  ylab("Plasma p-tau217/AB42")+
  scale_y_continuous(breaks = round(seq(0, 0.15, by = 0.05), 3), limits = c(0, 0.15)) +
  scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(0, 4),  labels = function(x) sprintf("%.2f", x)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"), guide = "none") + #003f5c
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 1.75, y = 0.125, label = "rho = 0.97, p = < 2.2e-16", size = 2)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Conc_ptau217_ab42.pdf",device = "pdf",width = 60, dpi=500, height = 55, units = "mm")

##p-tau217 and CSF

cor.test(df.val$ptau217, df.val$csf_4240, method="spearman")

ggplot(df.val, aes(x = ptau217, y = csf_4240, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.072, linetype="dotted", color="grey20", linewidth = 0.4)+
  geom_vline(xintercept=0.23211, linetype="dotted", color="grey20", linewidth = 0.4)+
  xlab("Plasma p-tau217 (pg/mL)")+
  ylab("CSF AB42/40")+
  scale_y_continuous(breaks = round(seq(0.0, 0.125, by = 0.05), 3), limits = c(0.0, 0.125)) +
  scale_x_continuous(breaks = seq(0, 4, by = 1), limits = c(0, 4),  labels = function(x) sprintf("%.2f", x)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"),guide="none")+
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 1.75, y = 0.100, label = "rho = -0.64, p = < 2.2e-16", size = 2)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Conc_ptau217_CSF.pdf",device = "pdf",width = 60, dpi=500, height = 55, units = "mm")

##p-tau217/AB42 and CSF

cor.test(df.val$PL_217_42, df.val$csf_4240, method="spearman")

ggplot(df.val, aes(x = PL_217_42, y = csf_4240, fill=petstat)) +
  geom_quasirandom(alpha = 0.9, size = 1, color="black", shape=21, stroke = 0.25)+
  geom_hline(yintercept=0.072, linetype="dotted", color="grey20", linewidth = 0.6)+
  geom_vline(xintercept=0.00738, linetype="dotted", color="grey20", linewidth = 0.6)+
  xlab("Plasma p-tau217/AB42")+
  ylab("CSF AB42/40")+
  scale_y_continuous(breaks = round(seq(0.0, 0.125, by = 0.05), 3), limits = c(0.0, 0.125)) +
  scale_x_continuous(breaks = round(seq(0, 0.15, by = 0.05), 3), limits = c(0, 0.15)) +
  scale_fill_manual(values = c("#c9c3c3", "#af639e"),guide="none")+
  theme_classic()+
  theme(axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  annotate("text", x = 0.1, y = 0.100, label = "rho = -0.66, p = < 2.2e-16", size = 2)
ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Conc_ptau217ab42_CSF.pdf",device = "pdf",width = 60, dpi=500, height = 55, units = "mm")



##in BFPC----

##in ADNI----

#S2 boxplots biomarkers and centiloids 24----
##BFMC----
df.val$petstat <- as.factor(ifelse(df.val$CL_fnc_ber_com_composite >=24, 1, 0))
df.val <- df.val %>% drop_na(petstat)
#("#d66f07", "#304C89","#ab123a")

ggplot(df.val, aes(x = petstat, y = ptau217, fill=petstat)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.08681, ymax = 0.23211,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  geom_hline(yintercept=0.23211, linetype="dotted", color="grey20") +
  geom_hline(yintercept=0.08681, linetype="dotted", color="grey20") +
  scale_y_continuous(breaks = round(seq(0, 2.2, by = 0.5), 2), limits = c(0, 2.2)) +
  scale_x_discrete(labels=c("0" ="<24 CL","1" ="≥24 CL")) +
  labs(x = "AD status", y = "Plasma p-tau217") +
  scale_fill_manual(values = c("#afafaf", "#ab123a"), guide = "none") +
  #scale_color_manual(values = c("#afafaf", "#ab123a"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 1)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Centiloids/LPptau217_Boxplot_CL24.pdf",device = "pdf",width = 50, dpi=500, height = 45, units = "mm")

ggplot(df.val, aes(x = petstat, y = PL_217_42, fill=petstat)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.0037, ymax = 0.0074,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  geom_hline(yintercept=0.00738, linetype="dotted", color="grey20") +
  geom_hline(yintercept=0.0037, linetype="dotted", color="grey20") +
  scale_y_continuous(breaks = round(seq(0, 0.085, by = 0.025), 3), limits = c(0, 0.085)) +
  scale_x_discrete(labels=c("0" ="<24 CL","1" ="≥24 CL")) +
  labs(x = "AD status", y = "Plasma p-tau217/AB42") +
  scale_fill_manual(values = c("#afafaf", "#304C89"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 0.04)
ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Centiloids/LPptau217AB42_Boxplot_CL24.pdf",device = "pdf",width = 50, dpi=500, height = 45, units = "mm")

ggplot(df.val, aes(x = petstat, y = csf_4240, fill=petstat)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.058, ymax = 0.073,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  geom_hline(yintercept=0.073, linetype="dotted", color="grey20") +
  geom_hline(yintercept=0.058, linetype="dotted", color="grey20") +
  scale_y_continuous(breaks = round(seq(0.025, 0.125, by = 0.025), 3), limits = c(0.025, 0.125)) +
  scale_x_discrete(labels=c("0" ="<24 CL","1" ="≥24 CL")) +
  labs(x = "AD status", y = "CSF AB42/40") +
  scale_fill_manual(values = c("#afafaf", "#d66f07"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold", size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 0.1)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/Centiloids/csf2_Boxplot_CL24.pdf",device = "pdf",width = 50, dpi=500, height = 45, units = "mm")


#S3 correlations between VR and PET----
##BFMC----
ggplot(df.val, aes(x = as.factor(VR_overall), y = CL_fnc_ber_com_composite, fill = as.factor(VR_overall)))+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.08681, ymax = 0.23211,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  geom_hline(yintercept=24, linetype="dotted", color="grey20") +
  #geom_hline(yintercept=0.08681, linetype="dotted", color="grey20") +
  scale_y_continuous(breaks = round(seq(-25, 175, by = 25), 1), limits = c(-25, 175)) +
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

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/Correlations_CL_VR.pdf",device = "pdf",width = 50, dpi=500, height = 45, units = "mm")

