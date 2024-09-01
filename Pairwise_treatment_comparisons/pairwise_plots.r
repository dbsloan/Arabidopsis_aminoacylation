library(ggplot2)

pairwise_data = read.csv("combined_pairwise_comparisons.csv")

pairwise_data$Genome = factor(pairwise_data$Genome, levels=c("Plastid", "Mitochondrial", "Nuclear"))

ggplot(data=pairwise_data, aes(x=CCA_1, y=CCA_2, color=Genome, size=0.1*Total_Reads)) + geom_abline(linewidth=0.25) + geom_point(alpha=0.3) + facet_grid (Comparison ~ Biorep) + theme_bw() + xlim(c(0,100)) + ylim(c(0,100)) + scale_color_brewer(palette = "Dark2") + theme(legend.position="top", legend.title=element_blank(), axis.text = element_text(size=7), axis.title = element_text(size=8), strip.text = element_text(size=8)) + xlab("Percentage CCA Tails") + ylab("Percentage CCA Tails") + scale_size(range = c(0.5, 3), guide="none")

#output: pairwise_treatment_correlation.pdf