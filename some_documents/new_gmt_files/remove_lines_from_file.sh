#!/bin/bash

i=156
while IFS=$'\t', read -r GSE Module_number
do
	if [ $GSE != "GSE" ]; then
		name="$GSE#$Module_number"
		echo $name
		index=$(grep -Pn "$name\t" "mm.new_$i.gmt" | grep -Po "\d+:")
		if [ ! -z "$index" ]; then
			length_index=$(echo -n "$index" | wc -c)
			index_in_gmt_file=${index:0:length_index-1}
			echo $index_in_gmt_file
			new_i=$(($i+1))
			echo "$index_in_gmt_filed"
			sed ${index_in_gmt_file}d "mm.new_$i.gmt">"mm.new_$new_i.gmt"
			i=$(($i+1))
		fi
	fi
done < ../new_modules_annotations_files/mm.wrong.tsv
