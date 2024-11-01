# End_mapping_CCA_vs_CC

Generating plots summarizing the start positions for reads from experiment 1 (periodate libraries only) for each isodecoder family, separating positions by reads that have intact CCA tail vs. those that end in CC.

`perl start_mapping_CC_CCA.pl file_list.all_reads.txt Arabidopsis_Sprinzl.mod.txt 100 > start_abundances.all_reads.txt`

Sprinzl coordinate file is found in "Sprinzl_coordinates" directory.

Output from above used to generate plots with split_CC_CCA_plots.all_reads.R. 