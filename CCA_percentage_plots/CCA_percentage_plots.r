library(ggplot2)
library (scales)

#Input uses total counts except for 5 pairs of organellar tRNAs with some ambiguous cross-mapping. Used unique counts for these

#mitochondrial-Ath-AsnGTT-183766	plastid-Ath-AsnGTT-171848
#mitochondrial-Ath-AspGTC-171828	plastid-Ath-AspGTC-171849
#mitochondrial-Ath-SerGGA-171839	plastid-Ath-SerGGA-171866
#mitochondrial-Ath-TrpCCA-171842	plastid-Ath-TrpCCA-171871
#mitochondrial-Ath-eMetCAT-171836	plastid-Ath-eMetCAT-171862

AA_data = read.csv("CCA.byAA.csv")

AA_data$Genome = factor(AA_data$Genome, levels=c("Plastid", "Mitochondrial", "Nuclear"))

#these values are from a global summing of all reads from all three libs. They are average CCA% by genome

hline_values = data.frame(Genome=c("Mitochondrial", "Nuclear", "Plastid"), Mean_CCA=c(62.9, 51.7, 62.2))

hline_values$Genome = factor(hline_values$Genome, levels=c("Plastid", "Mitochondrial", "Nuclear"))

ggplot(data=AA_data, aes(x=AA, y=100*CCA/(CCA+CC), color=Genome, shape=as.factor(BioRep))) + geom_point(alpha=0.5) + geom_hline(data=hline_values, aes(yintercept=Mean_CCA, color=Genome)) + theme_bw() + xlab ("tRNA Isoacceptor Families") + ylab ("Percentage CCA Tails") + ylim(c(0,100)) + facet_grid(~Genome, scales="free_x", space="free_x") + scale_color_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6.5), axis.text.y = element_text(size=7), legend.position = "none", axis.title=element_text(size=8), strip.text.x = element_text(size = 8))

#output: CCA.byAA.pdf


isodecoder_data=read.csv("CCA.byIsodecoder.csv")


ggplot(data=isodecoder_data, aes(x=Isodecoder, y=100*CCA/(CCA+CC), color=Genome, shape=as.factor(BioRep))) + geom_point(alpha=0.5, size=1) + geom_hline(data=hline_values, aes(yintercept=Mean_CCA, color=Genome)) + theme_bw() + xlab ("tRNA Isodecoder Families") + ylab ("Percentage CCA Tails") + ylim(c(0,100)) + facet_grid(~Genome, scales="free_x", space="free_x") + scale_color_brewer(palette = "Dark2") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5.25), axis.text.y = element_text(size=6), legend.position = "none", axis.title=element_text(size=8), strip.text.x = element_text(size = 8))

#output: CCA.byIsodecoder.pdf
