library(tidyverse)


CCA_data = read.csv("combined_pairwise_comparisons.aaRS_class.csv")

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
    ylim(c(0,1)) + xlab("Isoacceptor Family") +
    ylab("Ratio of CCA Percentages (Pre-deacylated : Periodate Treated)")

ggsave("~/Downloads/Predeacylation_AA_retention.pdf", width=7, height=2.5)