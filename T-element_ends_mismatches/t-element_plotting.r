library(ggplot2)

stemloop_combo = read.csv("stemloop_combo_for_plotting.csv")

stemloop_combo$Side = factor(stemloop_combo$Side, levels = c("Start", "End"))

ggplot(data=stemloop_combo, aes(x=Pos, y=Count, fill=Type)) + geom_hline(yintercept=0, color="gray80", linewidth=0.1) + geom_col() + facet_grid(Side ~ StemLoop, scales="free_y") + theme_bw() + ylab("Read Count") +xlab ("Read Mapping Position") + scale_fill_manual(values=c("gray20","firebrick")) + theme(legend.key.size = unit(0.3, "cm"), legend.position = "left", strip.text = element_text(size=6, margin = margin(2,2,2,2)), axis.text = element_text(size=6), legend.title = element_blank(), legend.text=element_text(size=5.5), axis.title = element_text(size=8,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggsave("t-element_end_pos.pdf", width=3.75, height=2.5)


mismatch_combo = read.csv("ccmC_t-element_IleCAT.csv")

mismatch_combo$MapStatus = factor(mismatch_combo$MapStatus, levels = c("Deletion", "A_Mismatch", "C_Mismatch", "G_Mismatch", "T_Mismatch", "Reference"))


ggplot(data=mismatch_combo, aes(x=Position, y=ReadCount, fill=MapStatus)) + geom_col() + facet_wrap(~Gene, nrow=2, scales="free_y")+ theme_bw() + scale_fill_manual(values = c("black","chartreuse4", "steelblue", "goldenrod", "firebrick", "gray70")) + ylab ("Read Count") + theme(legend.key.size = unit(0.3, "cm"), strip.text = element_text(size=6, margin = margin(2,2,2,2)), axis.text = element_text(size=6), legend.title = element_blank(), legend.text=element_text(size=5.5), axis.title = element_text(size=8,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("misincorp.pdf", height=2.25, width=5)
