# take base names of all files to combine all 4 lanes of sequencing into one file each read
while read -r basename; do
	cat ${basename}_*R1_* >> ../${basename}_R1.fastq.gz
	cat ${basename}_*R2_* >> ../${basename}_R2.fastq.gz
done < basenames.txt

# alternative method with numbered individuals
for i in {1..52}; do 
	for j in $( ls ${i}_*R1* ); do
 	cat $j >> /path/to/output/${i}_R1.fastq.gz; done
  	for k in $( ls ${i}_*R2* ); do
 	cat $k >> /path/to/output/${i}_R2.fastq.gz; done
done


