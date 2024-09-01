library(tidyverse)
library(RColorBrewer)

## set working directory and files there in
dir <- "."
filenames <- list.files(path=dir, pattern=c('CCA.txt'), full.names = TRUE)

## create dataframe for CCA vs non CCA trnas
table <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(table) <- c('Ref_tRNA', 'CCA_count', 'CC_count', "Other_count")

## upload Watkins Results
for (j in filenames){

  tablej <- read.csv(j, sep = "\t", header = T) %>% filter(Ref_tRNA != "Total")
  tablej$Lib <- j %>% str_replace(".*/", "")  %>% str_replace(".bbmerge.*", "")  
  
  table <- rbind(table,tablej)
  
}

table <- table %>% mutate(Genome = str_replace(Ref_tRNA, "_.*", ""))
table$Genome <- table$Genome %>% str_replace("Homo", "Nuclear") %>% str_replace("mito", "Mitochondrial")
table <- table %>% mutate(Nreads = CCA_count+CC_count+Other_count)
table$Isodecoders <- table$Ref_tRNA %>% str_replace("-[0-9]+.*", "") %>% str_replace(".*tRNA-", "")
table$Isoacceptor <- table$Isodecoders %>% str_replace("-[A-Z]+$", "") %>% str_replace(".*tRNA-", "") %>% str_replace("mito_", "") %>% str_replace("[1-2]+$", "")

## rename mt tRNAs
# Create a named vector for mapping full names to abbreviated names
tRNA_map <- c(
  TRNA = "Ala", TRNC = "Cys", TRND = "Asp", TRNE = "Glu", TRNF = "Phe", 
  TRNG = "Gly", TRNH = "His", TRNI = "Ile", TRNK = "Lys", TRNL = "Leu", 
  TRNM = "Met", TRNN = "Asn", TRNP = "Pro", TRNQ = "Gln", TRNR = "Arg", 
  TRNS = "Ser", TRNT = "Thr", TRNV = "Val", TRNW = "Trp", TRNY = "Tyr",
  Met = "Met(e)", iMet = "Met(i)"
)
table <- table %>%
  mutate(Isoacceptor = recode(Isoacceptor, !!!tRNA_map))


## Isoacceptors
tableIsoacceptor <- table %>% group_by(Isoacceptor, Lib, Genome) %>% summarise(CCA = sum(CCA_count), CC = sum(CC_count), Other = sum(Other_count))
tableIsoacceptor$Genome = factor(tableIsoacceptor$Genome, levels=c("Mitochondrial", "Nuclear"))

## Calculate average
Mean_CCAs <- tableIsoacceptor %>% group_by(Genome) %>% summarise(round(mean(100*CCA/(CCA+CC)), 1)) %>% select(2) %>% as.vector()
## Average horizontal lines
hline_values = data.frame(Genome=c("Mitochondrial", "Nuclear"), Mean_CCA= Mean_CCAs[[1]])
hline_values$Genome = factor(hline_values$Genome, levels=c("Mitochondrial", "Nuclear"))

# Define custom colors
custom_colors <- c("Mitochondrial" = scales::brewer_pal(palette = "Dark2")(3)[2],
                   "Nuclear" = scales::brewer_pal(palette = "Dark2")(3)[3])
## Plot
ggplot(data=tableIsoacceptor, aes(x=Isoacceptor, y=100*CCA/(CCA+CC), color=Genome, shape=as.factor(Lib))) + 
  geom_point(alpha=0.5) + 
  geom_hline(data=hline_values, aes(yintercept=Mean_CCA, color=Genome)) + 
  theme_bw() + 
  xlab ("tRNA Isoacceptor Families") + 
  ylab ("Percentage CCA Tails") + 
  ylim(c(0,100)) + 
  facet_grid(~Genome, scales="free_x", space="free_x") + 
  scale_color_manual(values = custom_colors) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6.5), axis.text.y = element_text(size=7), legend.position = "none", axis.title=element_text(size=8), strip.text.x = element_text(size = 8))
ggsave(paste(dir, "/CCA_tail_percentage.pdf", sep = ""), width = 4.5, height = 2.25)
