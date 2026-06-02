#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=test
#SBATCH --partition nocona
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-46

source activate bcftools

threads=2

# define main working directory
workdir=/lustre/scratch/jmanthey/10_ant_test

# base name of fastq files, intermediate files, and output files
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir}/basenames.txt | tail -n1 )

# define the reference genome
refgenome=/home/jmanthey/denovo_genomes/camp_sp_genome_filtered.fasta

# define the location of the reference mitogenomes
mito=/home/jmanthey/denovo_genomes/formicinae_mitogenomes.fasta

# define the location of the reference blochmannia genomes
bloch=/home/jmanthey/blochmannia/combined.fasta

# run bbduk
/lustre/work/jmanthey/bbmap/bbduk.sh \
in1=${workdir}/00_fastq/${basename_array}_R1.fastq.gz \
in2=${workdir}/00_fastq/${basename_array}_R2.fastq.gz \
out1=${workdir}/01_cleaned/${basename_array}_R1.fastq.gz \
out2=${workdir}/01_cleaned/${basename_array}_R2.fastq.gz \
minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 \
ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

source activate hist_filter

# seqprep adapter trimming and merging short inserts that were sequenced
SeqPrep -f ${workdir}/01_cleaned/${basename_array}_R1.fastq.gz \
-r ${workdir}/01_cleaned/${basename_array}_R2.fastq.gz \
-1 ${workdir}/01_cleaned/${basename_array}_unmerged_R1.fastq.gz \
-2 ${workdir}/01_cleaned/${basename_array}_unmerged_R2.fastq.gz \
-y J \
-s ${workdir}/01_cleaned/${basename_array}_merged.fastq.gz

source activate bcftools

# run bwa mem for unmerged reads
bwa-mem2 mem -t ${threads} ${refgenome} \
${workdir}/01_cleaned/${basename_array}_unmerged_R1.fastq.gz \
${workdir}/01_cleaned/${basename_array}_unmerged_R2.fastq.gz > \
${workdir}/01_bam_files/${basename_array}_unmerged.sam

# run bwa mem for merged reads
bwa-mem2 mem -t ${threads} ${refgenome} \
${workdir}/01_cleaned/${basename_array}_merged.fastq.gz > \
${workdir}/01_bam_files/${basename_array}_merged.sam

# convert sam to bam for unmerged reads
samtools view -b -S -@ ${threads} \
-o ${workdir}/01_bam_files/${basename_array}_unmerged.bam \
${workdir}/01_bam_files/${basename_array}_unmerged.sam

# convert sam to bam for merged reads
samtools view -b -S -@ ${threads} \
-o ${workdir}/01_bam_files/${basename_array}_merged.bam \
${workdir}/01_bam_files/${basename_array}_merged.sam

# combine the two bam files
samtools cat ${workdir}/01_bam_files/${basename_array}_merged.bam \
${workdir}/01_bam_files/${basename_array}_unmerged.bam > \
${workdir}/01_bam_files/${basename_array}.bam

# clean up the bam file
java -jar /home/jmanthey/picard.jar CleanSam \
I=${workdir}/01_bam_files/${basename_array}.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file
java -jar /home/jmanthey/picard.jar SortSam \
I=${workdir}/01_bam_files/${basename_array}_cleaned.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# add read groups to sorted and cleaned bam file
java -jar /home/jmanthey/picard.jar AddOrReplaceReadGroups \
I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam \
O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
RGLB=1 RGPL=illumina RGPU=unit1 RGSM=${basename_array}

# remove cleaned and sorted bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam

# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
java -jar /home/jmanthey/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
M=${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt \
I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam \
O=${workdir}/01_bam_files/${basename_array}_final.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted_rg.bam

# index the final bam file
samtools index ${workdir}/01_bam_files/${basename_array}_final.bam

# alignment stats
echo ${basename_array} > ${basename_array}.stats

# samtools depth sum of aligned sites
echo "samtools depth sum of aligned sites" >> ${basename_array}.stats
samtools depth  ${workdir}/01_bam_files/${basename_array}_final.bam  |  \
awk '{sum+=$3} END { print "Sum = ",sum}' >> ${basename_array}.stats

# proportion dupes
echo "proportion duplicates" >> ${basename_array}.stats
head -n8 ${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt | \
tail -n1 | cut -f9 >> ${basename_array}.stats

# samtools flagstat

samtools flagstat ${workdir}/01_bam_files/${basename_array}_final.bam > ${basename_array}.flagstat




