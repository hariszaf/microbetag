#!/bin/bash

# Input file
input_file="bins_taxa_original.csv"

# Folder containing files to be renamed
folder_path="bins"

# Loop through each file in the folder
for filename in "$folder_path"/*.fa; do

    # Extract the second column from the filename
    bin_name=$(basename "$filename" ".fa")
    echo ""
    echo "Bin name" $bin_name

    # Find the corresponding entry in the input file
    corresponding_bin=$(grep "$bin_name" "$input_file" | awk -F"\t" '{print $1}')

    if [ -n "$corresponding_bin" ]; then
	echo "===corresponding: " "$corresponding_bin"
	new_filename="$folder_path/${corresponding_bin}.fa"
	echo $new_filename
	mv "$filename" "$new_filename"
    else
        echo "!!! NO GOOD:" "$filename"
    fi
done

