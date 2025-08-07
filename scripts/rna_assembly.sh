minimap2 -ax splice -G 10k Necha2_scaffolds.fasta Necha2_ESTs.fasta | samtools sort > Necha2_aln.bam
stringtie -o Necha2_transcripts.gtf Necha2_aln.bam
