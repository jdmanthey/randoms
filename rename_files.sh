# tab delimited file with original names \t new names
while read -r name1 name2; do
	mv $name1 $name2
done < rename_text.txt
