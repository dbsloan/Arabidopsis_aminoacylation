library(tidyverse)
library(scales)

tpm_data = read.csv("t-element_vs_mito_tpms.csv")

ggplot(data=tpm_data, aes(x=Ref_tRNA, y=TPM,0.1, shape=as.factor(Rep), color=Type)) +
  geom_point(alpha=0.5, size=0.75) +
  facet_wrap(~Type, scales = "free_x") +
  theme_bw() +
  xlab("tRNA or tRNA-like Sequence") +
  ylab("Transcripts per Million") +
  theme(legend.position = "none", axis.title=element_text(size=7, face="bold"), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6), axis.text.y=element_text(size=6), strip.text = element_text(size = 7)) +
  scale_color_manual(values=c("black", "gray40")) +
  scale_y_log10(labels = label_number())

ggsave("t-element_vs_mito_tpms.pdf", width=3.25, height=2.25)