#!/usr/bin/env Rscript

library(tidyverse)
library(RColorBrewer)
library(reshape2)


## load table
edition <- read.csv("base-modification-table.txt", sep = "\t", header = T) 

## Plots
dir <- "."

##DBS: rewrote the original code (commented out below) to plot the total misincorp freq for each site, rather than a separate point for each type of mismatch per site.
total_freq = edition %>% rowwise() %>% mutate(Freq = sum(A_perc, G_perc, C_perc, T_perc, Del_perc, na.rm=TRUE))
total_freq <- total_freq %>% group_by(GeneName, Position, Reference, Compartment, Replicate) %>% filter(any(Freq>10))
total_freq <- dcast(total_freq, GeneName + Position + Reference + Compartment + Replicate ~ Tail, value.var = "Freq")
total_freq <- total_freq2 %>% filter(!is.na(CC) & !is.na(CCA))

total_freq$Compartment = factor(total_freq$Compartment, levels=c("plastid", "mitochondrial", "nuclear"), labels=c("Plastid", "Mitochondrial", "Nuclear"))

total_freq %>% ggplot(aes(x=CCA, y=CC, color = Reference)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~Compartment) +
  theme_bw() + 
  xlab ("Misincorporation Percentage in CCA Reads") + 
  ylab ("Misincorporation Percentage in CC Reads") + 
  xlim(c(0,100)) + 
  ylim(c(0,100)) + 
  scale_color_manual(values = c("chartreuse4", "steelblue", "goldenrod", "firebrick")) + 
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7), axis.title=element_text(size=8), strip.text.x = element_text(size = 8)) +
  guides(color = guide_legend(title = "Reference Base", label.position = "right"))

ggsave(paste(dir, "/CC-CCA_base-modification-periodate.pdf", sep = ""), width = 7.5, height = 2.5)


## base edition patterns in organellar genomes
#long_edition <- melt(edition, id.vars = c("GeneName", "Position", "Reference", "Compartment","Replicate",  "Tail", "Total_wDel"), measure.vars = c( "A_perc", "C_perc", "G_perc", "T_perc", "Del_perc"), variable.name = "base_mod", value.name = "Freq") 
#long_edition <- long_edition %>% filter(!is.na(Freq))
## remove bases with both CCA and CC < 10% edition
#long_edition <- long_edition %>% group_by(GeneName, Position, Reference, Compartment, Replicate) %>% filter(any(Freq>10))

## correlation of base-modification frequencies between aminoacylated and no-aminoacylated tRNAs
#long_edition <- dcast(long_edition, GeneName + Position + Reference + Compartment + Replicate + base_mod ~ Tail, value.var = "Freq")
#long_edition <- long_edition %>% filter(!is.na(CC) & !is.na(CCA))
#long_edition %>% ggplot(aes(x=CCA, y=CC, color = base_mod)) +
#  geom_point(alpha = 0.5) +
#  geom_abline(intercept = 0, slope = 1) +
#  facet_wrap(~Compartment) +
#  theme_bw() + 
#  xlab ("CCA") + 
#  ylab ("CC") + 
#  xlim(c(0,100)) + 
#  ylim(c(0,100)) + 
#  scale_color_brewer(palette = "Dark2") + 
#  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7), axis.title=element_text(size=8), strip.text.x = element_text(size = 8)) +
#  ggtitle("Base modification frequencies on Periodate Treatment", subtitle = "CC+CCA > 100 reads per pos // Mapping: nuclear = all, organellar = unique)") +
#  guides(color = guide_legend(title = "Modification", label.position = "right"))
#ggsave(paste(dir, "/CC-CCA_base-modification-periodate-no-thsd.pdf", sep = ""), units = "cm", width = 25, height = 12)
