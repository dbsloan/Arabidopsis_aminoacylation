library(tidyverse)

IleTAT = read.csv("nuclear-Ath-IleTAT.depth.combo.csv")

IleTAT_long = IleTAT %>% pivot_longer(cols = c("Reference", "A_Mismatch", "C_Mismatch", "G_Mismatch", "T_Mismatch", "Deletion"), names_to = "MapStatus", values_to = "ReadCount")

IleTAT_long$MapStatus = factor(IleTAT_long$MapStatus, levels = c("Deletion", "A_Mismatch", "C_Mismatch", "G_Mismatch", "T_Mismatch", "Reference"))

ggplot(data=IleTAT_long, aes(x=POS, y=ReadCount, fill=MapStatus)) + geom_col() + facet_wrap(~Trt, nrow=2, scales="free_y")+ theme_bw() + scale_fill_manual(values = c("black","chartreuse4", "steelblue", "goldenrod", "firebrick", "gray70")) + xlab ("Position") + ylab ("Read Count") + theme(legend.key.size = unit(0.3, "cm"), strip.text = element_text(size=6, margin = margin(2,2,2,2)), axis.text = element_text(size=6), legend.title = element_blank(), legend.text=element_text(size=5.5), axis.title = element_text(size=8,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("IleTAT.pdf", width=3.8, height=2.5)

LysTTT = read.csv("plastid-Ath-LysTTT.depth.combo.csv")

LysTTT_long = LysTTT %>% pivot_longer(cols = c("Reference", "A_Mismatch", "C_Mismatch", "G_Mismatch", "T_Mismatch", "Deletion"), names_to = "MapStatus", values_to = "ReadCount")

LysTTT_long$MapStatus = factor(LysTTT_long$MapStatus, levels = c("Deletion", "A_Mismatch", "C_Mismatch", "G_Mismatch", "T_Mismatch", "Reference"))

ggplot(data=LysTTT_long, aes(x=POS, y=ReadCount, fill=MapStatus)) + geom_col() + facet_wrap(~Trt, nrow=2, scales="free_y")+ theme_bw() + scale_fill_manual(values = c("black","chartreuse4", "steelblue", "goldenrod", "firebrick", "gray70")) + xlab ("Position") + ylab ("Read Count") + theme(legend.key.size = unit(0.3, "cm"), strip.text = element_text(size=6, margin = margin(2,2,2,2)), axis.text = element_text(size=6), legend.title = element_blank(), legend.text=element_text(size=5.5), axis.title = element_text(size=8,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("LysTTT.pdf", width=3.8, height=2.5)