#!/usr/bin/env Rscript

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)

#load Srinzl positions
sprinzl <- read.csv("../Sprinzl_coordinates/Arabidopsis_Sprinzl.mod.txt", sep = "\t", header = T) %>% select(-6)
sprinzl <- sprinzl %>% group_by(GeneName, Region) %>% mutate(id = row_number())

sprinzl <- sprinzl %>% mutate(Region = if_else(SprinzlPos == 12 | SprinzlPos == 13, "5' D stem", Region))
sprinzl <- sprinzl %>% mutate(Region = if_else(SprinzlPos == 14 | SprinzlPos == 21, "D loop", Region))
sprinzl <- sprinzl %>% mutate(Region = if_else(SprinzlPos == 22 | SprinzlPos == 23, "3' D stem", Region))


## Load Sprinzl positions vector
id2 <- sprinzl %>% filter(grepl("x", SprinzlPos)) %>% group_by(GeneName, Region) %>% mutate(id2 = row_number())
sprinzl <- left_join(sprinzl, id2)
#rename Sprinzl positions adding  suffixes numbering
sprinzl$SprinzlPos <- if_else(!is.na(sprinzl$id2), paste(sprinzl$SprinzlPos, sprinzl$id2, sep = ""), paste(sprinzl$SprinzlPos, sep = ""))
sprinzl <- sprinzl %>% select(-id2)

## Assign each tRNA position a letter 
RegionsOrder <- data.frame(
  Region = c("P -1", "5' acceptor stem", "P 8-9", "5' D stem", "D loop", "3' D stem", "P 26", "5' anticodon stem", "Anticodon loop", "3' anticodon stem", "Variable loop", "5' T stem", "T loop", "3' T stem", "3' acceptor stem", "Disciminator base", "CCA tail" ),
  letter = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q")
)

## add them to dataframe
sprinzl <- left_join(sprinzl, RegionsOrder)
## Join Region (letter) and Sprinzl position with suffix numbering
sprinzl <- sprinzl %>% mutate(letter = paste(letter, SprinzlPos, sep = " "))

    ## Function to sort tRNA regions including suffixes
    custom_sort <- function(pos_vec) {
      df <- data.frame(pos_vec = pos_vec) %>%
        # Separate the letter part and the numeric part (with suffixes)
        separate(pos_vec, into = c("letter_part", "numeric_suffix_part"), sep = " ") %>%
        mutate(
          numeric_part = as.integer(gsub("[^0-9].*$", "", numeric_suffix_part)),
          suffix_part = gsub("^[0-9]+", "", numeric_suffix_part),
          suffix_num_part = as.integer(gsub("^[^0-9]+", "", suffix_part)),
          suffix_num_part = ifelse(is.na(suffix_num_part), 0, suffix_num_part) # Handle cases without numeric suffix
        ) %>%
        # Sort by letter, numeric part, and numeric suffix part
        arrange(letter_part, numeric_part, suffix_num_part, suffix_part) %>%
        # Recombine the sorted parts into the original format
        mutate(sorted_pos_vec = paste(letter_part, numeric_suffix_part)) %>%
        pull(sorted_pos_vec)
      
      return(df)
    }



## load table
base_mod <- read.csv("tRNA-base-modification-average.txt", sep = "\t", header = T)

## add Sprinzl positions
base_mod <- left_join(base_mod, sprinzl)

#DBS: renaming Met genes and plastid/chloroplasts for more convenient alphabetical sorting 
base_mod$Isodecoder <- base_mod$Isodecoder %>% str_replace("eMetCAT", "Met(e)CAT")
base_mod$Isodecoder <- base_mod$Isodecoder %>% str_replace("iMetCAT", "Met(i)CAT")
base_mod$GeneName <- base_mod$GeneName %>% str_replace("eMetCAT", "Met(e)CAT")
base_mod$GeneName <- base_mod$GeneName %>% str_replace("iMetCAT", "Met(i)CAT")
base_mod$GeneName <- base_mod$GeneName %>% str_replace("plastid", "chloroplast")

## DBS: I am removing the non-canonical "x" Sprinzl position. They are almost never modified, most genes don't have them, and they were defined in a somewhat sloppy fashion in my original Sprinzl script. 
base_mod <- base_mod %>% filter(!str_detect(SprinzlPos, "x"))
sprinzl <- sprinzl %>% filter(!str_detect(SprinzlPos, "x"))

## sort Sprinzl positions
levels_plot <- sprinzl$letter %>% unique() %>% custom_sort()



## re-shape dataframe to have genes as rows and sprinzl positions as columns, filled by base modification percentages
data <- dcast(base_mod, GeneName + Position + SprinzlPos ~ letter, value.var = "base_mod_total_avg") 
data[is.na(data)] <- 0
## get only one rown per gene
collapsed_data <- data %>%
  group_by(GeneName) %>%
  summarise(across(3:(ncol(data)-1), sum, na.rm = TRUE))

## get all sprinzl positions as vector
new_columns <- sprinzl$letter %>% unique() %>% as.vector()
## add absent sprinzl positions
for (col in new_columns) {
  if (!col %in% colnames(collapsed_data)) {
    collapsed_data[[col]] <- 0
  }
}

## set colors for mt, nu, cp rows
gene_col <- collapsed_data %>% select(1) %>% unique()
gene_col$Genome <- gene_col$GeneName %>% str_replace("-.*", "")
gene_col<- gene_col %>% column_to_rownames('GeneName')


#set colors for columns
region_col <- collapsed_data %>% colnames() %>% as.data.frame()
colnames(region_col) <- "region"
region_col$letter <- region_col$region %>% str_replace(" .*", "")
region_col<- left_join(region_col, RegionsOrder) %>% select(-2)
region_col<- region_col %>% column_to_rownames('region')
region_col$Region <- region_col$Region %>% str_replace("3' ", "") %>% str_replace("5' ", "")

base_palette <- brewer.pal(11, "Spectral")
base_palette <- colorRampPalette(base_palette)(13)


annotation_colors <- list(
  Genome = c("nuclear" = scales::brewer_pal(palette = "Dark2")(3)[3],
             "chloroplast" = scales::brewer_pal(palette = "Dark2")(3)[1], 
             "mitochondrial" = scales::brewer_pal(palette = "Dark2")(3)[2]),
  Region = c("P -1" = base_palette[1], 
             "acceptor stem" = base_palette[2], 
             "P 8-9" = base_palette[3], 
             "D stem" = base_palette[4], 
             "D loop" = base_palette[5], 
             "P 26" = base_palette[6], 
             "anticodon stem" = base_palette[7], 
             "Anticodon loop" = base_palette[8], 
             "Variable loop" = base_palette[9], 
             "T stem" = base_palette[10], 
             "T loop" = base_palette[11], 
             "Disciminator base"= base_palette[12], 
             "CCA tail"= base_palette[13] )
)

## Prepare matrix to draw
## Rename rows and delete first column
collapsed_data <- collapsed_data %>% column_to_rownames('GeneName')
## sort columns in dataframe
collapsed_data <- collapsed_data[, levels_plot]
## transform dataframe to matrix
collapsed_data_matrix <- as.matrix(collapsed_data)

custom_colors <- colorRampPalette(c("ivory", "dodgerblue4"))(100)
pdf("heatmap.pdf", width = 15, height = 14)
pheatmap(
  collapsed_data_matrix,
  annotation_row = gene_col,
  annotation_colors = annotation_colors,
  annotation_col = region_col,
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  display_numbers = FALSE, 
  color = custom_colors, 
  fontsize_col = 8, 
  fontsize_row = 5,
  main = "Base modification frequency in control treatment"
)
dev.off()

