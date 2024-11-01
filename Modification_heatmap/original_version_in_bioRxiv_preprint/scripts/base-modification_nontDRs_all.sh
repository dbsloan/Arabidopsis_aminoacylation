#!/bin/sh

#SBATCH --time=24:00:00
#SBATCH --job-name=edition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12


for file in *read-mapping-positions.txt; do

	mkdir "${file%%.read-mapping-positions.txt}"-edition-nontDRs

	## put references in folder
	cp combined_Ath_tRNAs_withStemLoops.fas "${file%%.read-mapping-positions.txt}"-edition-nontDRs/combined_Ath_tRNAs_withStemLoops.fa

	#get nontDR read names
	awk '$5 != "tDR"' $file | cut -f 1 > "${file%%.read-mapping-positions.txt}"-edition-nontDRs/"${file%%.read-mapping-positions.txt}".nontDRs-reads
	
	cd "${file%%.read-mapping-positions.txt}"-edition-nontDRs
	
	cp ../"${file%%.read-mapping-positions.txt}"*blast_processed.fas blast-processed.tmp
	perl -i -pe 's/ .*//g' blast-processed.tmp

	#extract nontDRs reads
	seqkit grep -nf "${file%%.read-mapping-positions.txt}".nontDRs-reads blast-processed.tmp > "${file%%.read-mapping-positions.txt}".nontDRs-reads.fas
	#map them again separetly
	bowtie2-build combined_Ath_tRNAs_withStemLoops.fa combined_Ath_tRNAs_withStemLoops
	for file in *.fas; do
	    bowtie2 -a --no-unal -p 12 -D 25 -R 5 -N 1 -L 12 -f --end-to-end -x combined_Ath_tRNAs_withStemLoops -U $file -S "${file%%.fas}".tRNAs-SLs.sam > "${file%%.fas}".tRNAs-SLs.bowtie2.log.txt 2> "${file%%.fas}".tRNAs-SLs.bowtie2.err.txt 
	    samtools view --threads 12 -h -b "${file%%.fas}".tRNAs-SLs.sam | samtools sort > "${file%%.fas}".tRNAs-SLs.bam
	    samtools index "${file%%.fas}".tRNAs-SLs.bam
	    bam-readcount -w 1 -f combined_Ath_tRNAs_withStemLoops.fa "${file%%.fas}".tRNAs-SLs.bam > "${file%%.fas}".tRNAs-SLs.rc
	    parse_bamreadcount_v2.sh "${file%%.fas}".tRNAs-SLs.rc 
	done
	rm *bt2

	cd ..

done
