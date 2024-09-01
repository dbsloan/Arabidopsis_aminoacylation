# CCA_mapping_runs

Raw fastq were processed as follows

- R1 and R2 with bbmerge
- Flanking sequence from MSR-seq adapters trimmed (MSR-seq_trim.p)
- Collapsed identical sequences (collapse_identical_seqs.pl)
- Blast collapsed sequences against reference database and trim 5' end of individual reads (not collapsed) based on blast output (trim_5prime_blast.pl)
- Map individually blast-trimmed reads to reference database with bowtie2
- Summarize counts for both total mapping reads and uniquely mapping reads by end position (CCA, CC, or other) (CC_vs_CCA_counter.sam4.pl)

Reference databases (fasta format) and custom scripts are in reference_databases and scripts subdirectories, respectively. Output files and bash run scripts are in experiment subdirectories.