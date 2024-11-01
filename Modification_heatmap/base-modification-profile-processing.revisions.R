#!/usr/bin/env Rscript

library(tidyverse)
library(RColorBrewer)
library(reshape2)

#load tRNA names and lengths
trnalengths <- read.csv("reference-lengths.txt", sep = "\t", header = F) 
colnames(trnalengths) <- c("GeneName", "length")

#load Srinzl positions
sprinzl <- read.csv("../Sprinzl_coordinates/Arabidopsis_Sprinzl.mod.txt", sep = "\t", header = T) %>% select(-6)
sprinzl <- sprinzl %>% group_by(GeneName, Region) %>% mutate(id = row_number())


## UNIQUE MAPPINGS
## set working directory and CCA edition files there in
dir <- "unique"
filenames <- list.files(path= dir, pattern=c('parsed.txt'), full.names = TRUE)

## create dataframe for CCA  edition
data <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(data) <- c("Lib", "GeneName",	"Position",	"Reference", "A", "C", "G", "T", "Ins", "Del", "Total")

## upload all files with coverages of mappping against the genomes
for (j in filenames){

  dataj <- read.csv(j, sep = "\t", header = T)  %>% select(-11)
  dataj$Lib <- j %>% str_replace(".*/", "") %>% str_replace(".bb.*", "")
  colnames(dataj) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")
  data <- rbind(data,dataj)
}


data$Mapping <- "unique"

base_mod <- data

## ALL mappings
## Set working directory
dir <- "all"
filenames <- list.files(path= dir, pattern=c('parsed.txt'), full.names = TRUE)

## create dataframe for CCA  base_mod
data <- data.frame(matrix(ncol = 14, nrow = 0))
colnames(data) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")

## upload all files with coverages of mappping against the genomes
for (j in filenames){

  dataj <- read.csv(j, sep = "\t", header = T) %>% select(-11)
  dataj$Lib <- j %>% str_replace(".*/", "") %>% str_replace(".bb.*", "")
  colnames(dataj) <- c("GeneName",	"Position",	"Reference",	"A", "C", "G", "T", "Total", "Ins", "Del", "Lib")
  data <- rbind(data,dataj)
}


data$Mapping <- "all"


## Merge unique and all mapping data
base_mod <- rbind(base_mod, data)


base_mod$Isodecoder <- base_mod$GeneName %>% str_replace(".*Ath-", "") %>% str_replace("-.*", "")
base_mod$Treatment <- base_mod$Lib %>% str_replace(".*At.", "") %>% str_replace("a", "periodate") %>% str_replace("c", "no_periodate") %>% str_replace("b", "periodate_deacylated") 
base_mod$Compartment <- base_mod$GeneName %>% str_replace("-.*", "")
base_mod$Replicate <- base_mod$Lib %>% str_replace("MSR_At", "") %>% str_replace(".$", "")
base_mod$Total_wDel <- base_mod$Total + base_mod$Del


## Select unique mappings for organellar genomes and all for nuclear genome
base_mod <- base_mod %>% filter((Compartment == "nuclear" & Mapping == "all" ) | 
                              (Compartment != "nuclear" & Mapping == "unique"))

##DBS: exclude mitochondrial stemloops
base_mod <- base_mod %>% filter(Isodecoder != "stemloop")  
                                  

##Select no-periodate treatment and positions with >50 reads (including deletions)
##DBS: lowered cutoff to 50.
base_mod <- base_mod %>% 
  group_by(GeneName, Lib, Mapping, Position) %>%
  filter(sum(Total_wDel)>50) %>%
  filter(Treatment == "no_periodate") 


## Calculate percentages
columns <- c('A', 'G', 'C', 'T', 'Del')
for (col in columns) {
  percent_col <- paste0(col, "_perc")
  base_mod[[percent_col]] <- round((base_mod[[col]] / base_mod$Total_wDel) * 100, 2)
}



## Calculate total base modification percent
base_mod <- base_mod %>%
  mutate(base_mod_total = case_when(
    Reference == "A" ~ (100-A_perc),
    Reference == "G" ~ (100-G_perc),
    Reference == "T" ~ (100-T_perc),
    Reference == "C" ~ (100-C_perc),
    TRUE ~ NA_real_  # Default case, if none of the above match
  ))


## Delete non-edited/reference percentages
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


base_mod <- base_mod %>% select(11,14,16,12,15,13,1,2,3,4,5,6,7,10,17:23)

## count the number of Replicates that have a position editied in a gene
base_mod <- base_mod %>% group_by(GeneName, Position) %>% mutate( nRep = n())
base_mod_avg <- base_mod %>% select(6:9) %>% unique()
## Calculate base_mod averages across Replicates
base_mod_avg <- base_mod %>% 
  filter(nRep > 1) %>% 
  group_by(GeneName, Position) %>% 
  summarize(A_avg = round(mean(A_perc), 2), 
            G_avg = round(mean(G_perc), 2), 
            C_avg = round(mean(C_perc), 2),
            T_avg = round(mean(T_perc), 2), 
            Del_avg = round(mean(Del_perc), 2),
            base_mod_total_avg = round(mean(base_mod_total), 2)) %>% 
  right_join(base_mod_avg, .)

##LFC:Add sprinzl positions
base_mod_avg <- left_join(base_mod_avg, sprinzl)
base_mod_avg$Compartment <- base_mod_avg$GeneName %>% str_replace("-.*", "")

##LFC:merge base modification averages for isodecoder families considering Sprinzl positions
base_mod_avg_isod <- base_mod_avg %>% group_by(Isodecoder, Compartment, SprinzlPos) %>%
  summarise(
    base_mod_total_avg = round(mean(base_mod_total_avg), 2) 
  )

## Export tables
dir <- "."
write.table(base_mod, paste(dir, "/tRNA-base-modification.txt", sep = ""), quote = F, row.names = F, sep = "\t")
write.table(base_mod_avg, paste(dir, "/tRNA-base-modification-average.txt", sep = ""), quote = F, row.names = F, sep = "\t")
write.table(base_mod_avg_isod, paste(dir, "/tRNA-base-modification-average-isodecoders.txt", sep = ""), quote = F, row.names = F, sep = "\t")