#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool



baseCommand: [ collapse_table.py ]



inputs:

   input_table:
    type: string?  #trimmomatic-phred.yaml#phred?
    inputBinding:
      prefix: -phred
      separate: false
      position: 4
    label: 'quality score format'
    doc: >
      Either PHRED "33" or "64" specifies the base quality encoding. Default: 64   


# Use as in the following examples:
#        ./collapse_table.py -i taxonomic_otu_table.biom -g FAPROTAX_database.txt -o functional_otu_table.tsv -r report.txt --group_leftovers_as 'other' --normalize_collapsed 'columns_after_collapsing' -v
#        ./collapse_table.py -i taxonomic_otu_table.tsv --groups_file FAPROTAX_database.txt -f -o functional_otu_table.tsv -r report.txt --column_names_are_in last_comment_line  --keep_header_comments --non_numeric consolidate -v --row_names_are_in_column "taxonomy" --omit_columns 0 --normalize_collapsed columns_before_collapsing --group_leftovers_as 'other'






doc: >
   FAPROTAX is a manually constructed database that maps prokaryotic taxa (e.g. genera or species) to metabolic or other ecologically relevant functions (e.g. nitrification, denitrification or fermentation), based on the literature on cultured representatives. For example, if all cultured species within a bacterial genus (or more precisely, all type strains of species) have been identified as denitrifiers, FAPROTAX assumes that all uncultured members of that genus are also denitrifiers. Functions represented in FAPROTAX focus on marine and lake biogeochemistry, particularly sulfur, nitrogen, hydrogen and carbon cycling, although other functions (e.g. plant pathogeneicity) are also included. The complete list of functional groups covered by FAPROTAX, as well as all literature used, can be found within the database itself. FAPROTAX comes with a versatile script (collapse_table.py) for converting prokaryotic taxon abundance profiles ("OTU tables" or "taxon tables") into putative functional group abundance profiles ("function tables"), based on the taxa identified in a sample and their functional annotations in FAPROTAX.

label: Functional Annotation of Prokaryotic Taxa (FAPROTAX)

$namespaces:
 edam: http://edamontology.org/
 s: http://schema.org/

$schemas:
 - http://edamontology.org/EDAM_1.16.owl
 - https://schema.org/version/latest/schemaorg-current-http.rdf

s:license: "https://www.gnu.org/licenses/gpl-3.0.en.html"

s:copyrightHolder:
    - name: "Haris Zafeiropoulos"
    - url: "https://hariszaf.github.io/"