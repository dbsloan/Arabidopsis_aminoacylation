library(ggplot2)

ap_trizol_data = read.csv("AP_trizol_comparison.csv")

ap_trizol_data2 = ap_trizol_data %>% filter(Genome != "Spike-in control")

ap_trizol_data2$Genome = factor(ap_trizol_data2$Genome, levels=c("Plastid", "Mitochondrial", "Nuclear", "Spike-in control"))


ggplot(data=ap_trizol_data2, aes(x=100*Acid_Phenol, y=100*Trizol, color=Genome)) + geom_abline(linewidth=0.25) + geom_point(alpha=0.3, size=1) + facet_grid (Count_Type ~ Trt) + theme_bw() + xlim(c(0,100)) + ylim(c(0,100)) + scale_color_brewer(palette = "Dark2") + theme(legend.position="top", legend.title=element_blank(), axis.text = element_text(size=7), axis.title = element_text(size=8), strip.text = element_text(size=8)) + xlab("Percentage in Acid-Phenol Preps") + ylab("Percentage in Trizol Preps")

#exported as AP_trizol_comparison2.pdf


spike_in_data=read.csv("spike-in_AP_trizol_totals.csv")

ggplot(data=spike_in_data, aes(x=Trt, y=Percentage, color=tRNA, group=tRNA)) + geom_point(alpha=0.5, size=1, position=position_dodge(width=0.2)) + facet_grid(rows=vars(spike_in_data$Measure), cols=vars(spike_in_data$Extraction)) + theme_bw() + ylim(c(0,100)) + scale_color_manual(values=c("gray20","firebrick")) + xlab("") + theme(legend.position = "top", axis.text=element_text(size=7), strip.text=element_text(size=7), axis.title=element_text(size=8), legend.title=element_blank(), legend.text=element_text(size=7)) + scale_x_discrete(labels = label_wrap(10))


#output: spike-in_AP_trizol.pdf