#!/usr/bin/env bash

BASE="/home/haris/Desktop/programs/"

cd $BASE/kofam_scan/

for faa in $BASE/microbetag/app/ref-dbs/gtdb_genomes/genomes_faa/*; do 

    echo $faa ; 
    date;
    ./exec_annotation -f mapper -o $BASE/microbetag/app/ref-dbs/gtdb_genomes/annotations/"${faa##*/}".tsv $faa;
    mv $faa $BASE/microbetag/app/ref-dbs/gtdb_genomes/tmp
    date;
done





