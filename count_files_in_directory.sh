
for i in $( ls -d */ ); do
	echo $i;
	find $i -type f | wc -l;
done


