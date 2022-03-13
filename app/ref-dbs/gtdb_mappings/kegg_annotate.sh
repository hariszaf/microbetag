#!/usr/bin/env bash

BASE="/home/haris/Desktop/programs/"

cd $BASE/kofam_scan/

for faa in $BASE/microbetag/app/ref-dbs/gtdb_genomes/test/*; do 

    echo $faa ; date; 
    #echo "${faa##*/}"
    ./exec_annotation -f mapper -o $BASE/microbetag/app/ref-dbs/gtdb_genomes/annotations/"${faa##*/}".tsv $faa;

done





