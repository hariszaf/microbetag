#! /usr/bin/gawk -f

# USAGE:
# ./map_accession_to_ncbi_id.awk taxmap_embl_ssu_ref_132.txt access2species_silva132.tsv > species2ncbiId.tsv


BEGIN {
   -F"\t"
}

(ARGIND=1) {
   accessio_to_ncbi_id[$1"."$2"."$3] = $NF
}

(ARGIND=2)  {
   if ($1 in accessio_to_ncbi_id) {
      
      print substr($0, index($0,$2)) "\t" accessio_to_ncbi_id[$1]
   }
}
