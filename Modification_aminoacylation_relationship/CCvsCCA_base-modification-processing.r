#!/usr/bin/env Rscript

library(tidyverse)
library(RColorBrewer)
library(reshape2)

#load tRNA names and lengths
trnalengths <- read.csv("../Modification_heatmap/reference-lengths.txt", sep = "\t", header = F) 
colnames(trnalengths) <- c("GeneName", "length")

#load Srinzl positions
sprinzl <- read.csv("../Sprinzl_coordinates/Arabidopsis_Sprinzl.mod.txt", sep = "\t", header = T) %>% select(-6)
sprinzl <- sprinzl %>% group_by(GeneName, Region) %>% mutate(id = row_number())


## UNIQUE MAPPINGS
## set working directory and CCA base_mod files there in
dir <- "unique"
filenames <- list.files(path= dir, pattern=c('CCAreads'), full.names = TRUE)

## create dataframe for CCA  base_mod
CCAedit <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(CCAedit) <- c("Lib", "GeneName",	"Position",	"Reference", "A", "C", "G", "T", "Ins", "Del", "Total")

## upload all files with coverages of mappping against the genomes
for (j in filenames){

  CCAeditj <- read.csv(j, sep = "\t", header = T)  %>% select(-11)
  CCAeditj$Lib <- j %>% str_replace(".*/", "") %>% str_replace(".bb.*", "")
  colnames(CCAeditj) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")
  CCAedit <- rbind(CCAedit,CCAeditj)
}

CCAedit$Tail <- "CCA"

## load CC base_mod files there in
filenames <- list.files(path= dir, pattern=c('CCreads'), full.names = TRUE)

## create dataframe for CCA  base_mod
CCedit <- data.frame(matrix(ncol = 14, nrow = 0))
colnames(CCedit) <- c("Lib", "GeneName",	"Position",	"Reference", "A", "C", "G", "T", "Total", "Ins", "Del")

## upload all files with coverages of mappping against the genomes
for (j in filenames){

  CCeditj <- read.csv(j, sep = "\t", header = T)  %>% select(-11)
  CCeditj$Lib <- j %>% str_replace(".*/", "") %>% str_replace(".bb.*", "")
  colnames(CCeditj) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")
  CCedit <- rbind(CCedit,CCeditj)
}

CCedit$Tail <- "CC"


CCedit$Mapping <- "unique"
CCAedit$Mapping <- "unique"

base_mod <- rbind(CCedit, CCAedit) 

## ALL mappings
## Set working directory
dir <- "all"
filenames <- list.files(path= dir, pattern=c('CCAreads'), full.names = TRUE)

## create dataframe for CCA  base_mod
CCAedit <- data.frame(matrix(ncol = 14, nrow = 0))
colnames(CCAedit) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")

## upload all files with coverages of mappping against the genomes
for (j in filenames){

  CCAeditj <- read.csv(j, sep = "\t", header = T) %>% select(-11)
  CCAeditj$Lib <- j %>% str_replace(".*/", "") %>% str_replace(".bb.*", "")
  colnames(CCAeditj) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")
  CCAedit <- rbind(CCAedit,CCAeditj)
}

CCAedit$Tail <- "CCA"

## load CC base_mod files there in
filenames <- list.files(path= dir, pattern=c('CCreads'), full.names = TRUE)

## create dataframe for CCA  base_mod
CCedit <- data.frame(matrix(ncol = 14, nrow = 0))
colnames(CCedit) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")

## upload all files with coverages of mappping against the genomes
for (j in filenames){

  CCeditj <- read.csv(j, sep = "\t", header = T) %>% select(-11)
  CCeditj$Lib <- j %>% str_replace(".*/", "") %>% str_replace(".bb.*", "")
  colnames(CCeditj) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")
  CCedit <- rbind(CCedit,CCeditj)
}

CCedit$Tail <- "CC"

CCedit$Mapping <- "all"
CCAedit$Mapping <- "all"


## Merge unique and all mapping data
base_mod <- rbind(base_mod,CCedit, CCAedit)

base_mod$Isodecoder <- base_mod$GeneName %>% str_replace(".*Ath-", "") %>% str_replace("-.*", "")
base_mod$Treatment <- base_mod$Lib %>% str_replace(".*At.", "") %>% str_replace("a", "periodate") %>% str_replace("c", "no_periodate") %>% str_replace("b", "periodate_deacylated") 
base_mod$Compartment <- base_mod$GeneName %>% str_replace("-.*", "")
base_mod$Replicate <- base_mod$Lib %>% str_replace("MSR_At", "") %>% str_replace(".$", "")
base_mod$Total_wDel <- base_mod$Total + base_mod$Del


## Select unique mappings for organellar genomes and all for nuclear genome
base_mod <- base_mod %>% filter((Compartment == "nuclear" & Mapping == "all" ) | 
                              (Compartment != "nuclear" & Mapping == "unique"))


##Select priodate treatment
base_mod <- base_mod %>% 
  group_by(GeneName, Lib, Tail, Mapping, Position) %>%
  filter(sum(Total)>100) %>%
  filter(Treatment == "periodate") 


## Calculate percentages
columns <- c('A', 'G', 'C', 'T', 'Del')
for (col in columns) {
  percent_col <- paste0(col, "_perc")
  base_mod[[percent_col]] <- round((base_mod[[col]] / base_mod$Total_wDel) * 100, 2)
}

## Define a minimum base modification threshold
#DBS: note that this does not actually filter because the threshhold value is defined as a string in quotes. Not numerical
thsd <- "10"
base_mod <- base_mod %>%
  filter(
    (Reference == "A" & (`C_perc` > thsd | `G_perc` > thsd | `T_perc` > thsd)) |
    (Reference == "C" & (`A_perc` > thsd | `G_perc` > thsd | `T_perc` > thsd)) |
    (Reference == "G" & (`A_perc` > thsd | `C_perc` > thsd | `T_perc` > thsd)) |
    (Reference == "T" & (`A_perc` > thsd | `C_perc` > thsd | `G_perc` > thsd))
  )

## Delete non-modified/reference percentages
base_mod <- base_mod %>%
  mutate(
    `A_perc` = case_when(
      Reference == "A" ~ NA_real_,
      TRUE ~ `A_perc`
    ),
    `C_perc` = case_when(
      Reference == "C" ~ NA_real_,
      TRUE ~ `C_perc`
    ),
    `G_perc` = case_when(
      Reference == "G" ~ NA_real_,
      TRUE ~ `G_perc`
    ),
    `T_perc` = case_when(
      Reference == "T" ~ NA_real_,
      TRUE ~ `T_perc`
    )
  )


base_mod <- base_mod %>% select(12,11,15,17,13,16,14,1,2,3,4,5,6,7,10,18:23)

## Export table
dir <- "."
write.table(base_mod, paste(dir, "/base-modification-table.txt", sep = ""), quote = F, row.names = F, sep = "\t")
