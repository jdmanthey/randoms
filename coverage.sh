# samtools depth mean of aligned sites
samtools depth  *bamfile*  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'

# samtools depth sum of aligned sites
samtools depth  *bamfile*  |  awk '{sum+=$3} END { print "Sum = ",sum}'

# samtools depth mean of all sites
samtools depth  -a *bamfile*  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
