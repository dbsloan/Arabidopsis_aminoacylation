library(tidyverse)
library(RColorBrewer)

## upload Arabidopsis data
tableArab <- read.csv("Total_counts.csv", sep = ",", header = T)

## set working directory and files there in
dir <- "."
filenames <- list.files(path=dir, pattern=c('CCA.txt'), full.names = TRUE)

## create dataframe for CCA vs non CCA trnas
table <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(table) <- c('Ref_tRNA', 'CCA_count', 'CC_count', "Other_count")

## upload Watkins Data
for (j in filenames){

  tablej <- read.csv(j, sep = "\t", header = T) %>% filter(Ref_tRNA == "Total")
  tablej$Lib <- j %>% str_replace(".*/", "")  %>% str_replace(".bbmerge.*", "")  
  
  table <- rbind(table,tablej)
  
}


table <- table %>% select(-1)
table$Trt <- "Watkins" 
table <-  mutate( table, CCA_Percentage = round(CCA_count*100 / (CCA_count + CC_count),1))
table <-  mutate( table, Other_Percentage = round(Other_count*100 / (CCA_count + CC_count+Other_count),1))

## Merge Arabidopsis and Watkins data
tableArab <- tableArab %>% select(-3)
colnames(table) <- colnames(table) %>% str_replace("_count", "")
table <- rbind(tableArab,table)


## Plot
ggplot(data=table, aes(x=Trt, y=CCA_Percentage)) + 
  geom_jitter(width=0.1, alpha=0.5) + 
  theme_bw() + 
  xlab ("Treatment") + 
  ylab ("Percentage CCA Tails") + 
  ylim(c(0,100))
ggsave(paste(dir, "/CCA_percentage.pdf", sep = ""), width = 10, height = 10, units = "cm")

ggplot(data=table, aes(x=Trt, y=Other_Percentage)) + 
  geom_jitter(width=0.1, alpha=0.5) + 
  theme_bw() + 
  xlab ("Treatment") + 
  ylab ("Non-CC/CCA Percentage") + 
  ylim(c(0,100))
ggsave(paste(dir, "/Other_percentage.pdf", sep = ""), width = 10, height = 10, units = "cm")
