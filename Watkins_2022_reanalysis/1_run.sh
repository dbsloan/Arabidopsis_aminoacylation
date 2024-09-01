#!/bin/sh

#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --error=slurm.stderr
#SBATCH --output=slurm.stdout

#Merge R1 and R2 seqs
for file in *1.fq.gz; do bbmerge.sh in1=$file \
	in2=${file%_1.fq.gz}_2.fq.gz \
	out=${file%_1.fq.gz}.bbmerge.fq \
	outu=${file%_1.fq.gz}.unmerged1.fq \
	outu2=${file%_1.fq.gz}.unmerged2.fq \
	ihist=${file%_1.fq.gz}.bbmerge_hist.txt \
	ordered=t \
	qtrim=r \
	minoverlap=30 \
	mismatches=4 2> ${file%_1.fq.gz}.bbmerge_log.txt
done

#trim adapter sequences
for file in SRR18302277.bbmerge.fq; do perl MSR-seq_trim.pl $file ACTACCA 6 > ${file%fq}trim.fas; done
for file in SRR18302278.bbmerge.fq; do perl MSR-seq_trim.pl $file ACTCAGA 6 > ${file%fq}trim.fas; done
for file in SRR18302279.bbmerge.fq; do perl MSR-seq_trim.pl $file ACTGGAA 6 > ${file%fq}trim.fas; done

#create fasta file with collapsed identical seqs
for file in SRR18302277*trim.fas; do perl collapse_identical_seqs.pl $file > ${file%trim.fas}collapsed.fas; done

#blast collapsed file against tRNA database and use blast output to generate a version of the (uncollapsed) read file excluding reads without hits and trimming 5' ends that run off the tRNA. The map reads with bowtie and parse sam file
makeblastdb -in combined_Hs_tRNAs_CCA.fa -dbtype nucl -parse_seqids
bowtie2-build combined_Hs_tRNAs_CCA.fa combined_Hs_tRNAs_CCA.fa

for file in *.bbmerge.collapsed.fas; do 
	blastn -task blastn -db combined_Hs_tRNAs_CCA.fa  -query $file -evalue 1e-12 -num_threads 12 -out ${file%collapsed.fas}blast.txt
	perl trim_5prime_blast.pl ${file%collapsed.fas}blast.txt $file ${file%collapsed.fas}trim.fas > ${file%collapsed.fas}blast_processed.fas
	bowtie2 --no-unal -p 12 -L 10 -i C,1 --mp 5,2 --score-min L,-0.7,-0.7 -f -x combined_Hs_tRNAs_CCA.fa -U ${file%collapsed.fas}blast_processed.fas -S ${file%collapsed.fas}sam
	perl CC_vs_CCA_counter.sam4.pl ${file%collapsed.fas}sam > ${file%collapsed.fas}CC_vs_CCA.txt
	perl CC_vs_CCA_counter.sam4.pl ${file%collapsed.fas}sam --unique > ${file%collapsed.fas}CC_vs_CCA.unique.txt
done
