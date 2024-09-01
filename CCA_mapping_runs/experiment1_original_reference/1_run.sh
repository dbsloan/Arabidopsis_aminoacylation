#!/bin/sh

#SBATCH --time=100000:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error=slurm.stderr
#SBATCH --output=slurm.stdout

#Merge R1 and R2 seqs
for file in *1.fq.gz; do bbmerge.sh in1=$file in2=${file%1.fq.gz}2.fq.gz out=${file%_CKDL230041724-1A_HHVLKDSX7_L4_1.fq.gz}.bbmerge.fq outu=${file%_CKDL230041724-1A_HHVLKDSX7_L4_1.fq.gz}.unmerged1.fq outu2=${file%_CKDL230041724-1A_HHVLKDSX7_L4_1.fq.gz}.unmerged2.fq ihist=${file%_CKDL230041724-1A_HHVLKDSX7_L4_1.fq.gz}.bbmerge_hist.txt ordered=t qtrim=r minoverlap=30 mismatches=4 2> ${file%_CKDL230041724-1A_HHVLKDSX7_L4_1.fq.gz}.bbmerge_log.txt; done

#trim adapter sequences
for file in *bbmerge.fq; do perl MSR-seq_trim.pl $file ACTGGAA 6 > ${file%fq}trim.fas; done

#create fasta file with collapsed identical seqs
for file in *trim.fas; do perl collapse_identical_seqs.pl $file > ${file%trim.fas}collapsed.fas; done

#blast collapsed file against tRNA database and use blast output to generate a version of the (uncollapsed) read file excluding reads without hits and trimming 5' ends that run off the tRNA. The map reads with bowtie and parse sam file
for file in *collapsed.fas; do blastn -task blastn -db combined_Ath_tRNAs_withStemLoopsAndBacillus.fas -query $file -evalue 1e-12 -num_threads 12 -out ${file%collapsed.fas}blast.txt; perl trim_5prime_blast.pl ${file%collapsed.fas}blast.txt $file ${file%collapsed.fas}trim.fas > ${file%collapsed.fas}blast_processed.fas; bowtie2 --no-unal -p 12 -L 10 -i C,1 --mp 5,2 --score-min L,-0.7,-0.7 -f -x combined_Ath_tRNAs_withStemLoopsAndBacillus.fas -U ${file%collapsed.fas}blast_processed.fas -S ${file%collapsed.fas}sam; perl CC_vs_CCA_counter.sam4.pl ${file%collapsed.fas}sam > ${file%collapsed.fas}CC_vs_CCA.txt; perl CC_vs_CCA_counter.sam4.pl ${file%collapsed.fas}sam --unique > ${file%collapsed.fas}CC_vs_CCA.unique.txt; done
