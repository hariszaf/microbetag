#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

# -i finalTable.tsv       --- input_table
# -o func_table.tsv       --- out_collapsed
# -g FAPROTAX.txt         --- input_groups_file
# -d "Classification"     --- row_names_are_in_column
# -c "#"                  --- comment_prefix
# -v                      --- verbose // does not need a value, it's a flag
# -s /mnt/                --- out_sub_tables_dir

requirements:
 ResourceRequirement:
   ramMin: 100240
   coresMin: 2
 ShellCommandRequirement: {}


baseCommand: [ collapse_table.py ]

inputs:

   input_table:
      type: File  
      format: edam:data_3707
      inputBinding:
         position: 1
      label: 'Path to OTU table to be collapsed'
      doc: >
         Path to a classical table file (TSV, CSV or similar) or a BIOM file. Row names are stored in a column specified by --row_names_are_in_column, the rest of the columns contain numerical or string data. Non-numerical entries will be interpreted according to --non_numeric. Alternatively, a .biom observation table can be given, in which case observation IDs will be taken as 'row names'.


   input_groups_file:
      type: File
      format: edam:format_2330
      inputBinding:
         position: 2
      label: 'Path to the FAPROTAX database file'
      doc: >
         Path to a file defining the groups by which to collapse the table. In the framework of microbetag this is always the FAPROTAX database .txt file
   

   row_names_are_in_column:
      type: String
      position: 4
      doc: >
          Column containing row names (as referred to in the groups unless --collapse_columns_instead_of_rows is set). If column names are available, this specifies a column by name, otherwise it specifies a column by index (starting at 0). If this is empty ('', default), then row names are set to row indices (starting at 0) and all data in the table is taken as is. Only relevant for classical input tables.

outputs:

   out_collapsed:
      type: File
      position: 3
      format: edam:format_3475
      label: 'Main FAPROTAX output file'
      outputBinding:
         glob: $(inputs.input_table.nameroot).funct_otu_table.tsv


arguments: 

   - prefix: -c
      valueFrom: '#'
      position: 
   - prefix: --out_sub_tables_dir
      valueFrom: '/mnt/'
      position: 
   - prefix: -v
      position

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