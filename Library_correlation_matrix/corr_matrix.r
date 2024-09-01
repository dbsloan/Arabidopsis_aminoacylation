library(tidyverse)
library(plotly)
library(GGally)

sum_isodecoders <- function(file_path) {
  df <- read.table(file_path, header = TRUE, sep = "\t")
  
  df_grouped = df %>% separate(col="Ref_tRNA", into = c("Genome", "Species", "Isodecoder", "GeneID"), sep = "-") %>% 
    filter(Genome != "Bacillus" & Genome != "Total" & Isodecoder != "stemloop") %>%
    group_by(Genome, Isodecoder) %>%
    summarize(
      CCA_count_sum = sum(CCA_count),
      CC_count_sum = sum(CC_count),
      Other_count_sum = sum(Other_count)
    )
  
  df_grouped
  
  sums <- rowSums(df_grouped[, 3:5])
  
  # Return the sum vector
  return(log10(sums))
}


Periodate1 = sum_isodecoders ("MSR_At1a.bbmerge.CC_vs_CCA.txt")
Periodate2 = sum_isodecoders ("MSR_At2a.bbmerge.CC_vs_CCA.txt")
Periodate3 = sum_isodecoders ("MSR_At3a.bbmerge.CC_vs_CCA.txt")
Predeacylated1 = sum_isodecoders ("MSR_At1b.bbmerge.CC_vs_CCA.txt")
Predeacylated2 = sum_isodecoders ("MSR_At2b.bbmerge.CC_vs_CCA.txt")
Predeacylated3 = sum_isodecoders ("MSR_At3b.bbmerge.CC_vs_CCA.txt")
NoPeriodate1 = sum_isodecoders ("MSR_At1c.bbmerge.CC_vs_CCA.txt")
NoPeriodate2 = sum_isodecoders ("MSR_At2c.bbmerge.CC_vs_CCA.txt")
NoPeriodate3 = sum_isodecoders ("MSR_At3c.bbmerge.CC_vs_CCA.txt")

log10_count_df = data.frame(NoPeriodate1, NoPeriodate2, NoPeriodate3, Periodate1, Periodate2, Periodate3, Predeacylated1, Predeacylated2, Predeacylated3)

ggpairs(log10_count_df, lower = list(continuous = wrap("points", size = 0.25, alpha=2))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), strip.text = element_text(size=5), axis.text = element_text(size=5), axis.title = element_text(size=7), strip.background = element_rect(size = 0.25)) + 
  xlab("log10 Read Count") + 
  ylab("log10 Read Count")

ggsave("corr_matrix_libraries/corr_matrix_plot.pdf", width=6.5, height=6.5)