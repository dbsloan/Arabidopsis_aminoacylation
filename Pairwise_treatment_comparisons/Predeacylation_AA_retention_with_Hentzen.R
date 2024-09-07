library(tidyverse)


CCA_data = read.csv("~/Documents/ColoradoState/projects/tRNAs/MSR-seq/Sequencing/20231206_novogene_At_MSR-seq/Analysis/20240620/aaRS_class/combined_pairwise_comparisons.aaRS_class.csv")

CCA_data = separate(CCA_data, col="Gene", into = c("Genome2", "Species", "AA", "ID_num"), sep = "-")
CCA_data$AA = substr(CCA_data$AA, 1, nchar(CCA_data$AA) - 3)
CCA_data$AA = CCA_data$AA %>% str_replace("eMet", "Met(e)")
CCA_data$AA = CCA_data$AA %>% str_replace("iMet", "Met(i)")
CCA_data $Genome = factor(CCA_data $Genome, levels=c("Plastid", "Mitochondrial", "Nuclear"))

CCA_data %>% filter(Comparison == "Periodate vs Deacylated") %>% 
  group_by(Genome, AA, Isodecoder) %>% 
  summarize(CCAratio = mean(CCA_2/CCA_1), Reads = mean(Total_Reads)) %>% 
  ggplot(aes(x=AA, y=CCAratio, color=Genome, size=Reads)) +
    geom_point(alpha=0.5) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_color_brewer(palette = "Dark2") +
    scale_y_log10(limits=c(0.01,1)) +
    xlab("Isoacceptor Family") +
    ylab("Ratio of CCA Percentages (Pre-deacylated : Periodate Treated)") +
    theme (axis.title = element_text(size=7), axis.text = element_text(size = 6), legend.position = "None")


ggsave("~/Downloads/Predeacylation_AA_retention_panelA.pdf", width=4.75, height=2.65)


half_life_data = read.csv("~/Documents/ColoradoState/projects/tRNAs/MSR-seq/Sequencing/20231206_novogene_At_MSR-seq/Analysis/20240620/aaRS_class/Hentzen1972.csv")

ggplot(half_life_data, aes (x=T.half, y=CCA_retention_ratio, label=AA)) +
  geom_smooth(method="lm", se=FALSE, size=0.5) +
  geom_point(alpha=0.5, shape=21, fill="black", stroke=NA) +
  geom_text(size=1.75, nudge_y=0.05) +
  theme_bw() +
  scale_y_log10(limits=c(0.01,1)) +
  scale_x_log10() +
  xlab ("Aminoacyl-tRNA Half Life (min) from Hentzen et al. 1972") +
  ylab ("CCA Retention Ratio") +
  theme (axis.title = element_text(size=7), axis.text = element_text(size = 6))

ggsave("~/Downloads/Predeacylation_AA_retention_panelB.pdf", width=2.25, height=2.5)


summary(lm (log10(CCA_retention_ratio) ~ log10(T.half), data=half_life_data))

cor.test(log10(half_life_data$CCA_retention_ratio), log10(half_life_data$T.half))
cor.test(half_life_data$CCA_retention_ratio, half_life_data$T.half)
