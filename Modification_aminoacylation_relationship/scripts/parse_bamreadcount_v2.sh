#!/bin/bash

bam_readcount=$1

awk -v OFS='\t' '{print $1,$2,$3,$6,$7,$8,$9,($6+$7+$8+$9)}' $bam_readcount |  
perl -pe 's/A\:(\d*).*\tC\:(\d*).*\tG\:(\d*).*\tT\:(\d*).*/$1\t$2\t$3\t$4/g' | awk 'BEGIN{FS=OFS="\t"} {print $0, ($4+$5+$6+$7)}' > "${bam_readcount%%.*}"_parsed.tmp 

echo -e "name\tpos\tref\tA\tC\tG\tT\tTotal\tIns\tDel" > head.tmp

cat head.tmp "${bam_readcount%%.*}"_parsed.tmp > "${bam_readcount%%.*}"_parsed.txt

#rm *.tmp




#extraigo las columnas con info de los indels
cat $bam_readcount | cut -f11- > indels.tmp
#cuento los indels por cada linea
python parse_bamindels-v2.py indels.tmp indels_parsed.tmp
#agrego las columnas de indels
paste "${bam_readcount%%.*}"_parsed.txt indels_parsed.tmp > "${bam_readcount%%.tRNAs-SLs.*}"_indels_parsed.txt


rm *.tmp
