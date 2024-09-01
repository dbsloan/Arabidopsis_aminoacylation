#!/bin/sh

#SBATCH --job-name=positions
#SBATCH --time=100000:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4


for file in *.sam; do
	samtools view -h --threads 4 $file | grep "^@SQ" | perl -pe 's/.*SN://g' | perl -pe 's/LN://g' > reference-lengths.txt
    sam2pos-v3.py $file --unique "${file%%.sam}".read-mapping-positions-unique.txt
    left-right-positions-count-tailed.py reference-lengths.txt "${file%%.sam}".read-mapping-positions-unique.txt "${file%%.sam}".map-pos-count-unique.txt
	sam2pos-v3.py $file "${file%%.sam}".read-mapping-positions.txt
	left-right-positions-count-tailed.py reference-lengths.txt "${file%%.sam}".read-mapping-positions.txt "${file%%.sam}".map-pos-count.txt
done
