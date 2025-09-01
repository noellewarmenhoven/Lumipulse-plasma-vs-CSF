#Figures manuscript

#1. Load data

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

df$petstat <- as.factor(df$VR_overall)

df.val <- df %>% filter(df$Study %in% "VALIDATE")
df.val <- df.val %>% select(sid_bf2, sid_diag, VR_overall, CL_fnc_ber_com_composite, csf_4240, ptau217, PL_217_42, Plasma_AB142_pgmL_Lumi_c, Plasma_AB142_pgmL_Lumi, PL_217_42_n, petstat)
df.val <- df.val %>% drop_na(csf_4240, ptau217)



#Fig 1: Boxplots----
ggplot(df.val, aes(x = petstat, y = ptau217, fill=petstat)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.19, ymax = 0.32,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  geom_hline(yintercept=0.19, linetype="dashed", color="grey20", alpha = .8,linewidth = .35) +
  geom_hline(yintercept=0.32, linetype="dashed", color="grey20", alpha = .8,linewidth = .35) +
  geom_hline(yintercept=0.253, linetype="dashed", color="red", alpha=.8,linewidth = .35) +
  scale_y_continuous(breaks = seq(0, 2.2, by = 0.50), limits = c(0, 2.2),
                     labels = function(x) sprintf("%.2f", x)) +
  scale_x_discrete(labels=c("0" ="VR-","1" ="VR+")) +
  labs(x = "AD status", y = "Plasma p-tau217 (pg/mL)") +
  scale_fill_manual(values = c("#afafaf", "#4ca481"), guide = "none") +
  #scale_color_manual(values = c("#afafaf", "#ab123a"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 1)

ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/BioFINDER/Validate/FDA approach/LPptau217_Boxplot.pdf",device = "pdf",width = 50, dpi=500, height = 45, units = "mm")

ggplot(df, aes(x = petstat, y = PL_217_42, fill=petstat)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.0037, ymax = 0.0074,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_hline(yintercept=0.00738, linetype="dashed", color="grey20", alpha=.8, linewidth = .35) +
  geom_hline(yintercept=0.0037, linetype="dashed", color="grey20", alpha = .8, linewidth = .35) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  scale_y_continuous(breaks = round(seq(0, 0.1, by = 0.025), 3), limits = c(0, 0.085)) +
  scale_x_discrete(labels=c("0" ="VR-","1" ="VR+")) +
  labs(x = "AD status", y = "Plasma p-tau217/AB42") +
  scale_fill_manual(values = c("#afafaf", "#c52851"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 0.04)
ggsave("~/Documents/Projects/LP vs CSF/LP vs CSF/AAIC/Validate/LPptau217AB42_Boxplot_PB_csf42ptau181.pdf",device = "pdf",width = 50, dpi=500, height = 45, units = "mm")


ggplot(df, aes(x = petstat, y = csf_4240, fill=petstat)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.058, ymax = 0.072,
           fill = "grey80", alpha = 0.4) +
  geom_violin(alpha=0.3) +
  geom_hline(yintercept=0.058, linetype="dashed", color="grey20", alpha=.8, linewidth = .35) +
  geom_hline(yintercept=0.072, linetype="dashed", color="grey20", alpha = .8, linewidth = .35) +
  geom_quasirandom(alpha = 0.75, size = 1,shape=21, color="black", stroke=0.25) +
  #scale_y_continuous(breaks = round(seq(0, 0.1, by = 0.025), 3), limits = c(0, 0.085)) +
  scale_x_discrete(labels=c("0" ="VR-","1" ="VR+")) +
  labs(x = "AD status", y = "CSF AB42/40") +
  scale_fill_manual(values = c("#afafaf", "#be5391"), guide = "none") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=7),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6))+
  stat_compare_means(method="t.test", label="p.signif", label.x = 1.5, label.y = 0.04)
