# take base names of all files to combine all 4 lanes of sequencing into one file each read
while read -r basename; do
	cat ${basename}*R1_* >> ../${basename}_R1.fastq.gz
	cat ${basename}*R2_* >> ../${basename}_R2.fastq.gz
done < basenames.txt

