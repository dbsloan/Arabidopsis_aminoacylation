#!/bin/sh

#SBATCH --time=100000:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error=slurm.b.stderr
#SBATCH --output=slurm.b.stdout

#input files were previously processed with bbmerge and upstream scripts (see experiment1_original_reference run)

for file in *b.bbmerge.collapsed.fas; do blastn -task blastn -db combined_Ath_tRNAs_withStemLoopsAndBacillus2.fas -query $file -evalue 1e-12 -num_threads 12 -out ${file%collapsed.fas}blast.txt; perl trim_5prime_blast.pl ${file%collapsed.fas}blast.txt $file ${file%collapsed.fas}trim.fas > ${file%collapsed.fas}blast_processed.fas; bowtie2 --no-unal -p 12 -L 10 -i C,1 --mp 5,2 --score-min L,-0.7,-0.7 -f -x combined_Ath_tRNAs_withStemLoopsAndBacillus2.fas -U ${file%collapsed.fas}blast_processed.fas -S ${file%collapsed.fas}sam; perl CC_vs_CCA_counter.sam4.pl ${file%collapsed.fas}sam > ${file%collapsed.fas}CC_vs_CCA.txt; perl CC_vs_CCA_counter.sam4.pl ${file%collapsed.fas}sam --unique > ${file%collapsed.fas}CC_vs_CCA.unique.txt; done
