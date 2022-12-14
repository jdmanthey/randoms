
# move all fastq to current directory
find . -name '*fastq.gz' -print0 | xargs -0 -I '{}' mv "{}" .
