library(tidyverse)

start_pos_data = read.table("start_abundances.all_reads.txt", header=TRUE)

start_pos_data$Isodecoder = start_pos_data$Isodecoder %>% str_replace("eMet", "Met(e)")
start_pos_data$Isodecoder = start_pos_data$Isodecoder %>% str_replace("iMet", "Met(i)")

start_pos_data$EndType = factor(start_pos_data$EndType, labels = c("CC (not aminoacylated)", "CCA (aminoacylated)"))

start_pos_data %>%
  filter(Genome == "nuclear") %>%
  ggplot(aes(x=Position, y=Abundance/3, fill=EndType)) +
    geom_col() +
    facet_wrap(~Isodecoder, ncol=6) +
    theme_bw() +
    scale_fill_manual(values = c("gray60", "black")) +
    xlim (c(-1,50)) +
    theme(legend.position = "top", legend.title = element_blank(), strip.text = element_text(size=7, margin=margin(t=1, b=1)), axis.title = element_text(size=7, face="bold"), axis.text = element_text(size=6), legend.text = element_text(size=7)) +
    ylab ("Mean Proportion of Reads")

ggsave("nuclear_start_CC_CCA_plot.all_reads.pdf", width=6.7, height=7.7)

start_pos_data %>%
  filter(Genome == "plastid") %>%
  ggplot(aes(x=Position, y=Abundance/3, fill=EndType)) +
  geom_col() +
  facet_wrap(~Isodecoder, ncol=6) +
  theme_bw() +
  scale_fill_manual(values = c("gray60", "black")) +
  xlim (c(-1,50)) +
  theme(legend.position = "top", legend.title = element_blank(), strip.text = element_text(size=7, margin=margin(t=1, b=1)), axis.title = element_text(size=7, face="bold"), axis.text = element_text(size=6), legend.text = element_text(size=7)) +
  ylab ("Mean Proportion of Reads")

ggsave("plastid_start_CC_CCA_plot.all_reads.pdf", width=6.7, height=5.2)

start_pos_data %>%
  filter(Genome == "mitochondrial") %>%
  filter(!Isodecoder %in% c("AspGTC", "SerGGA", "TrpCCA")) %>%
  ggplot(aes(x=Position, y=Abundance/3, fill=EndType)) +
  geom_col() +
  facet_wrap(~Isodecoder, ncol=6) +
  theme_bw() +
  scale_fill_manual(values = c("gray60", "black")) +
  xlim (c(-1,50)) +
  theme(legend.position = "top", legend.title = element_blank(), strip.text = element_text(size=7, margin=margin(t=1, b=1)), axis.title = element_text(size=7, face="bold"), axis.text = element_text(size=6), legend.text = element_text(size=7)) +
  ylab ("Mean Proportion of Reads")


ggsave("mitochondrial_start_CC_CCA_plot.all_reads.pdf", width=6.7, height=3.5)


