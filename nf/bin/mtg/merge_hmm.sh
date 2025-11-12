#!/usr/bin/env bash
threads="$1"
output_file="$2"

echo -e "bin_id\tcontig_id\tko_term" > $output_file

ls hmmout_*/*.hmmout \
  | parallel -j $threads --bar --no-notice '
      f={}

      # Skip file if it has no non-comment lines
      if ! grep -q "^[^#]" "$f"; then
          exit 0
      fi

      fname=$(basename "$f")                  # e.g., K00099_bin_48.hmmout
      dir=$(dirname "$f")                     # e.g., hmmout_bin_48
      dirbase=$(basename "$dir")              # e.g., hmmout_bin_48
      genome=${dirbase#hmmout_}                  # e.g., bin_48
      ko=${fname%%_*}                          # e.g., K00099

      # Process non-comment lines
      grep "^[^#]" "$f" | awk -F" " -v genome="$genome" -v ko="$ko" '\''{print genome "\t" $1 "\t" ko}'\''
  ' >> $output_file
