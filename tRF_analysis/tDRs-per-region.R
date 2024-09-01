#!/usr/bin/env Rscript

library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(Polychrome)

setwd(".")

#load tRNA names and lengths
trnalengths <- read.csv("../Modification_heatmap/reference-lengths.txt", sep = "\t", header = F) 
colnames(trnalengths) <- c("GeneName", "length")

#load Srinzl positions
sprinzl <- read.csv("../Sprinzl_coordinates/Arabidopsis_Sprinzl.mod.txt", sep = "\t", header = T) %>% select(-6)
sprinzl <- sprinzl %>% group_by(GeneName, Region) %>% mutate(id = row_number())

## set working directory of UNIQUE MAPPINGS and files there in
dir <- "unique"
filenames <- list.files(path= dir, pattern=c('unique.txt'), full.names = TRUE)

## create dataframe
poscount <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(poscount) <- c("Lib", "GeneName", "Position", "Left_Other", "Left_tDR", "Left_C", "Left_CC", "Left_CCA", "Right_Other", "Right_tDR", "Right_C", "Right_CC", "Right_CCA", "Left_Total", "Right_Total", "Mapping")

## upload all files with mapping position counts
for (j in filenames){
  poscountj <- read.csv(j, sep = "\t", header = T)
  #colnames(CCAtmp) <- c("Ref_tRNA","CCA_count", "CC_count", "C_count", "Other_count")
  poscountj$Lib <- j %>% str_replace(".*/", "") %>% str_replace(".bb.*", "")
  poscountj$Mapping <- "unique"
  poscount <- rbind(poscount,poscountj)
}


## set working directory of ALL MAPPINGS and files there in
dir <- "all"
filenames <- list.files(path= dir, pattern=c('.txt'), full.names = TRUE)

## upload all files with mapping position counts
for (j in filenames){
  poscountj <- read.csv(j, sep = "\t", header = T)
  #colnames(CCAtmp) <- c("Ref_tRNA","CCA_count", "CC_count", "C_count", "Other_count")
  poscountj$Lib <- j %>% str_replace(".*/", "") %>% str_replace(".bb.*", "")
  poscountj$Mapping <- "all"
  poscount <- rbind(poscount,poscountj)
}

## define new variables
poscount$Isodecoder <- poscount$GeneName %>% str_replace(".*Ath-", "") %>% str_replace("-.*", "")
poscount$Treatment <- poscount$Lib %>% str_replace(".*At.", "") %>% str_replace("a", "periodate") %>% str_replace("c", "no_periodate") %>% str_replace("b", "periodate_deacylated") 
poscount$Genome <- poscount$GeneName %>% str_replace("-.*", "")

##filter
poscount <- poscount %>% filter(!grepl("stemloop", GeneName))
poscount <- poscount %>% filter(Treatment == "no_periodate")

## filter All mappings for nuclear genes, and Unique mappings for organellar genes
poscount <- poscount %>% filter((Genome == "nuclear" & Mapping == "all" ) | 
                              (Genome != "nuclear" & Mapping == "unique"))

## getting totals reads per tRNA per library (sum(Left_Total) or sum(Right_Total) should be the same)
poscount <- poscount %>% group_by(GeneName, Lib) %>% mutate(TotalReads = sum(Left_Total)) %>% ungroup()
## filter tRNAs with >150 reads mapped
poscount <- poscount %>% filter(TotalReads>150)
  
## calculate number of tDRs per tRNA per Library
poscount <- poscount %>% group_by(GeneName, Lib) %>% mutate(tDRreads = sum(Right_tDR)) %>% ungroup()

## Calculate tDR frequency for each position compared to the number of reads mapped to tRNA
poscount$tDR_Freq <- round(poscount$Right_tDR/poscount$TotalReads, 3)

## add tRNA lengths
poscount <- left_join(poscount, trnalengths)

## add tRNA Sprinzl positions
poscount <- left_join(poscount, sprinzl, by = c("GeneName", "Position"))

## add position in variable nucleotides/loops to Sprinzl position 
id2 <- poscount %>% filter(grepl("x", poscount$SprinzlPos)) %>% group_by(GeneName, Region, Lib) %>% mutate(id2 = row_number())
id2 <- id2 %>% select(Lib,GeneName,Position,SprinzlPos,Nucleotide,Region,id, id2)
poscount <- left_join(poscount, id2)
poscount$SprinzlPos <- if_else(!is.na(poscount$id2), paste(poscount$SprinzlPos, poscount$id2, sep = ""), paste(poscount$SprinzlPos, sep = ""))
poscount <- poscount %>% select(-id2)

## Assign each tRNA position a letter 
RegionsOrder <- data.frame(
  Region = c("P -1", "5' acceptor stem", "P 8-9", "5' D stem", "D loop", "3' D stem", "P 26", "5' anticodon stem", "Anticodon loop", "3' anticodon stem", "Variable loop", "5' T stem", "T loop", "3' T stem", "3' acceptor stem", "Disciminator base", "CCA tail" ),
  letter = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q")
)

poscount <- left_join(poscount, RegionsOrder)
poscount <- poscount %>% mutate(letter = paste(letter, Region, sep = " "))

## Some positions with freq>0.05 are consecutives, which I understand as "wider peaks" due to degradation. Also, sometimes they start at one tRNA region and continues in the next.
## This function selects only the 3'-most position of each peak. Can consider peaks for which the 3'-most position are >distance_threshold far from each other (I used 3 bp)
select_highest_consecutive <- function(df, group_cols, value_col, distance_threshold) {
  group_syms <- syms(group_cols) # Convert group_cols to symbols
  
  df %>%
    group_by(across(all_of(group_cols))) %>%             # Group by the specified group columns
    arrange(across(all_of(group_cols)), !!sym(value_col)) %>% # Ensure values are in order within groups
    mutate(consecutive = c(FALSE, diff(!!sym(value_col)) == 1)) %>% # Identify consecutive values
    group_by(across(all_of(group_cols)), rleid = cumsum(!consecutive)) %>% # Group by groups and consecutive runs
    filter(if(any(consecutive)) !!sym(value_col) == max(!!sym(value_col)) else TRUE) %>% # Select the highest value among consecutives or keep all
    ungroup() %>%
    group_by(across(all_of(group_cols))) %>%             # Re-group by the specified group columns
    mutate(
      dist_from_prev = c(Inf, diff(!!sym(value_col))),   # Calculate distances from previous positions
      group_num = cumsum(dist_from_prev >= distance_threshold)
    ) %>%
    group_by(across(all_of(group_cols)), group_num) %>%
    filter(if(n() == 1) TRUE else !!sym(value_col) == max(!!sym(value_col))) %>%
    ungroup() %>%
    select(-consecutive, -rleid, -dist_from_prev, -group_num)  # Clean up temporary columns
}

poscount_peak_3most_position <- poscount %>% filter(tDR_Freq>0.05)
poscount_peak_3most_position <- select_highest_consecutive(poscount_peak_3most_position, c("GeneName", "Lib"), "Position", 3)

## PLOTS
## Frequency of tDRs
dir <- "."

poscount$Genome = factor(poscount$Genome, levels = c("plastid", "mitochondrial", "nuclear"), labels = c("Plastid", "Mitochondrial", "Nuclear"))

poscount %>% ggplot(aes(x = 100*tDRreads/TotalReads)) +
  geom_histogram(aes(y = 100 * ..density..), fill = "dodgerblue4" , bins = 30, alpha = 0.7) +
  geom_vline(xintercept = 100 * 0.05, color = "black", linetype = "dashed", size = 0.5) +
  theme_bw() +
  facet_grid(~Genome) +
  scale_x_continuous(breaks = seq(0,100, 20)) +
  xlab("tRF Read Percentage") +
  ylab("Frequency Distribution (%)") +
  theme(axis.text.x = element_text(size=7), 
        axis.text.y = element_text(size=7), 
        legend.position = "none", 
        axis.title=element_text(size=8), 
        strip.text.x = element_text(size = 8),
        panel.grid.minor = element_blank())

#  ggtitle("Frequency distribution of tDRs per tRNA", subtitle = "dashed line: minimum threshold to consider tDRs // reads per tRNA >150")
ggsave(paste(dir, "/tDR-proportion.pdf", sep = ""), device = "pdf", width = 6.7, height = 2.25)

## 3' end of tDR per tRNA region (isodecoders)
tDR_stats <- poscount_peak_3most_position %>% 
  group_by(GeneName, Isodecoder, Genome, SprinzlPos, letter) %>% 
  summarize(nReplicates = n()) %>% ungroup() %>% 
  filter(nReplicates >1)  
tDR_stats_count <- tDR_stats %>% group_by(Isodecoder,Genome, SprinzlPos, letter) %>% summarize(isodecoder_count = n())

tDR_stats_count$Genome = factor(tDR_stats_count$Genome, levels = c("plastid", "mitochondrial", "nuclear"), labels = c("Plastid", "Mitochondrial", "Nuclear"))

tDR_stats_count %>% ggplot(aes(x = letter, fill = letter)) +
  geom_bar() +
  theme_bw() +
  facet_grid(~Genome) +
  scale_fill_brewer(palette = "Dark2") + 
  ylab("Number of tRFs") +
  xlab("tRNA region") +
  scale_y_continuous(breaks = seq(0,20,5)) +
#  ggtitle("tDR ends (3')", subtitle = "filters: control, freq >0.05, reads per tRNA >150, replicates >1") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6.5), 
        axis.text.y = element_text(size=7), 
        legend.position = "none", 
        axis.title=element_text(size=8), 
        strip.text.x = element_text(size = 8))
ggsave(paste(dir, "/tDR-3-end-regions_isodecoder.pdf", sep = ""), device = "pdf", width = 6.7, height = 3)

## prepare table to export
colnames(tDR_stats) <- c("Gene Name", "Isodecoder", "Genome", "Sprinzl Position", "tRNA Region", "Replicates")
tDR_stats <- tDR_stats %>% select(3,2,1,4,5,6)
tDR_stats$`tRNA Region` <- tDR_stats$`tRNA Region`%>% str_replace("^[A-Z] ", "")
## export table
write.table(tDR_stats, paste(dir, "/tDR-3-end-regions_isodecoder_table.txt", sep = ""), sep = "\t", quote = F, row.names = F)
