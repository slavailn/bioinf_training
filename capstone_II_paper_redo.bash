# Search GEO for a bioproject and download runInfo
esearch -db sra -query PRJNA639351 | efetch -format runinfo

# Download the reads in SRA format
~/programs/sratoolkit.3.0.1-ubuntu64/bin/prefetch -v -p SRR12010076
~/programs/sratoolkit.3.0.1-ubuntu64/bin/prefetch -v -p SRR12010075
~/programs/sratoolkit.3.0.1-ubuntu64/bin/prefetch -v -p SRR12010077

# Convert .sra to .fastq
 ~/programs/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump --split-3 SRR12010075/SRR12010075.sra
 ~/programs/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump --split-3 SRR12010076/SRR12010076.sra
 ~/programs/sratoolkit.3.0.1-ubuntu64/bin/fasterq-dump --split-3 SRR12010077/SRR12010077.sra
 
# Initial quality control
fastqc *.fastq

# Download the reference sequence
esearch -db nucleotide -query "CP052388.1" | efetch -format fasta > Klebsiella_pneumoniae.fasta
esearch -db nucleotide -query "CP052388.1" | efetch -format genbank > Klebsiella_pneumoniae.gbk
# Install this tool to convert genebank to gff2
conda install -c bioconda readseq
readseq Klebsiella_pneumoniae.gbk -f 24 # gff is format 24, consult "format" section in help file

# Create simulated data, take a look at fastqc results for the original reads
# to get a rough idea about the error rates
dwgsim -N 1300000 -1 150 -2 150 -e 0.0001-0.002 -E 0.0001-0.002 ../genome/Klebsiella_pneumoniae.fasta output
# Check the quality of simulated data
fastqc output.bwa.read*.fastq.gz

# Map the reads to the reference genome using BWA
bwa index Klebsiella_pneumoniae.fasta
bwa mem -t 2 genome/bwa_index/Klebsiella_pneumoniae.fasta fastq/SRR12010075_1.fastq  fastq/SRR12010075_2.fastq > SRR12010075.sam
bwa mem -t 2 genome/bwa_index/Klebsiella_pneumoniae.fasta fastq/SRR12010076_1.fastq  fastq/SRR12010076_2.fastq > SRR12010076.sam
bwa mem -t 2 genome/bwa_index/Klebsiella_pneumoniae.fasta fastq/SRR12010077_1.fastq  fastq/SRR12010077_2.fastq > SRR12010077.sam
bwa mem -t 2 ../genome/bwa_index/Klebsiella_pneumoniae.fasta output.bwa.read1.fastq.gz  output.bwa.read2.fastq.gz > simulated.sam

# Convert to bam, sort, deduplicate and index
samtools view -@ 2 -S -b SRR12010075.sam > SRR12010075.bam
samtools view -@ 2 -S -b SRR12010076.sam > SRR12010076.bam
samtools view -@ 2 -S -b SRR12010077.sam > SRR12010077.bam
samtools view -@ 2 -S -b simulated.sam > simulated.bam

# Sort
samtools sort -@ 2 SRR12010075.bam -o SRR12010075.sorted.bam
samtools sort -@ 2 SRR12010076.bam -o SRR12010076.sorted.bam
samtools sort -@ 2 SRR12010077.bam -o SRR12010077.sorted.bam
samtools sort -@ 2 simulated.bam -o simulated.sorted.bam

# Markduplicates
sambamba markdup -t 2 SRR12010075.sorted.bam SRR12010075.dedup.bam
sambamba markdup -t 2 SRR12010076.sorted.bam SRR12010076.dedup.bam
sambamba markdup -t 2 SRR12010077.sorted.bam SRR12010077.dedup.bam
sambamba markdup -t 2 simulated.sorted.bam simulated.dedup.bam

# Index
samtools index SRR12010075.dedup.bam
samtools index SRR12010076.dedup.bam
samtools index SRR12010077.dedup.bam
sambamba index simulated.dedup.bam

# Check alignment qualities
qualimap bamqc -nt 2 -bam align/SRR12010075.dedup.bam -outdir bamqc_results -outfile SRR12010075.html
qualimap bamqc -nt 2 -bam align/SRR12010076.dedup.bam -outdir bamqc_results -outfile SRR12010076.html
qualimap bamqc -nt 2 -bam align/SRR12010077.dedup.bam -outdir bamqc_results -outfile SRR12010077.html
qualimap bamqc -nt 2 -bam simulated.dedup.bam -outdir bamqc_results -outfile simulated.html

# Call variants with bcftools
bcftools mpileup --threads 2 -Ovu -f genome/Klebsiella_pneumoniae.fasta align/SRR12010075.dedup.bam > SRR12010075.vcf
bcftools mpileup --threads 2 -Ovu -f genome/Klebsiella_pneumoniae.fasta align/SRR12010076.dedup.bam > SRR12010076.vcf
bcftools mpileup --threads 2 -Ovu -f genome/Klebsiella_pneumoniae.fasta align/SRR12010077.dedup.bam > SRR12010077.vcf
bcftools mpileup --threads 2 -Ovu -f ../genome/Klebsiella_pneumoniae.fasta simulated.dedup.bam > simulated.gt.vcf

bcftools call --ploidy 1 -vm -O v SRR12010075.vcf > SRR12010075.call.vcf
bcftools call --ploidy 1 -vm -O v SRR12010076.vcf > SRR12010076.call.vcf
bcftools call --ploidy 1 -vm -O v SRR12010077.vcf > SRR12010077.call.vcf
bcftools call --ploidy 1 -vm -O v simulated.gt.vcf > simulated.call.vcf

# Collect BCF stats
bcftools stats SRR12010075.call.vcf > SRR12010075.call.stats
bcftools stats SRR12010075.call.vcf > SRR12010075.call.stats
bcftools stats SRR12010075.call.vcf > SRR12010075.call.stats
bcftools stats simulated.call.vcf > simulated.call.stats

# Visualize vcf files and aligments in IGV

# Compare real dataseq to ground truth
parallel -j 2 "bgzip {}" ::: *.call.vcf # bgzip vcf files
parallel -j 2 "bcftools index {}" ::: *.call.vcf.gz
parallel -j 2 "tabix {}" ::: *.call.vcf.gz

bgzip simulated.call.vcf
bcftools index simulated.call.vcf.gz
tabix simulated.call.vcf.gz

bgzip output.mutations.vcf
bcftools index output.mutations.vcf.gz
tabix output.mutations.vcf.gz

# Run a comparison: simulated data to ground truth
bcftools isec simulated.call.vcf.gz output.mutations.vcf.gz -p vcf_comparison
vcf-compare simulated.call.vcf.gz  output.mutations.vcf.gz > sim_vs_truth.txt

# Try a different approach with freebayes
freebayes -f genome/Klebsiella_pneumoniae.fasta -p 1 sim_data/simulated.dedup.bam > freebayes.vcf
bcftools index freebayes.vcf.gz
tabix freebayes.vcf.gz
vcf-compare freebayes_calls/freebayes.vcf.gz  output.mutations.vcf.gz > freebayes_vs_truth.txt # compare to ground truth

# Filter the variants in a similar fashion as in Tagliani et al. 2017:
# Read depth >= 30; Allele Freq >= 0.5; supporting forwar reads > 4 and supporting reverese reads > 4
SnpSift filter "(DP >= 30) & (AF >= 0.5) & (SAR > 4) & (SAF > 4)" freebayes_calls/freebayes.vcf.gz > freebayes.tagliani.vcf
bgzip freebayes.tagliani.vcf
tabix freebayes.tagliani.vcf.gz
vcf-compare freebayes.tagliani.vcf.gz  output.mutations.vcf.gz

# BONUS TASK: call the variants using snippy and compare to ground truth
https://github.com/tseemann/snippy







