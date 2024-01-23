

awk -F"\t"  'NR==FNR { map[$2]=$1; next } { if ($1 in map) $1=map[$1]; if ($2 in map) $2=map[$2]; print $1, $2, $3 }' gtdb2patricIds.tsv overlaps_07.tsv > gtdb_overlaps_07.tsv 



