#!/usr/bin/env cwl-runner

## REMEMBER!!! This works if and only if you are running cwl-runner from /home/haris 
## Make the appropriate edits to move it on microbetag's Docker environment. 

cwlVersion: v1.2
class: CommandLineTool
label: Functional Annotation of Prokaryotic Taxa (FAPROTAX)
doc: |
  FAPROTAX is a manually constructed database that maps prokaryotic taxa (e.g. genera or species) to metabolic or other ecologically relevant functions (e.g. nitrification, denitrification or fermentation), based on the literature on cultured representatives. For example, if all cultured species within a bacterial genus (or more precisely, all type strains of species) have been identified as denitrifiers, FAPROTAX assumes that all uncultured members of that genus are also denitrifiers. Functions represented in FAPROTAX focus on marine and lake biogeochemistry, particularly sulfur, nitrogen, hydrogen and carbon cycling, although other functions (e.g. plant pathogeneicity) are also included. The complete list of functional groups covered by FAPROTAX, as well as all literature used, can be found within the database itself. FAPROTAX comes with a versatile script (collapse_table.py) for converting prokaryotic taxon abundance profiles ("OTU tables" or "taxon tables") into putative functional group abundance profiles ("function tables"), based on the taxa identified in a sample and their functional annotations in FAPROTAX.
$namespaces:
  edam: http://edamontology.org/
  s: http://schema.org/

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: 2
    ramMin: 100240

inputs:
  input_groups_file:
    type: File
    default:
      class: File
      location: ./FAPROTAX_1.2.4/FAPROTAX.txt
    inputBinding:
      prefix: --input_groups_file
      position: 4
  input_table:
    label: Path to OTU table to be collapsed
    doc: |
      Path to a classical table file (TSV, CSV or similar) or a BIOM file. Row names are stored in a column specified by --row_names_are_in_column, the rest of the columns contain numerical or string data. Non-numerical entries will be interpreted according to --non_numeric. Alternatively, a .biom observation table can be given, in which case observation IDs will be taken as 'row names'.
    type: File
    format: edam:data_3707
    inputBinding:
      prefix: --input_table
      position: 2
  output_directory_path:
    type: string
    default: faprotax_output
  row_names_are_in_column:
    doc: |
      Column containing row names (as referred to in the groups unless --collapse_columns_instead_of_rows is set). If column names are available, this specifies a column by name, otherwise it specifies a column by index (starting at 0). If this is empty ('', default), then row names are set to row indices (starting at 0) and all data in the table is taken as is. Only relevant for classical input tables.
    type: string
    inputBinding:
      prefix: --row_names_are_in_column
      position: 5
  scrit:
    type: File
    default:
      class: File
      location: ./collapse_table.py
    inputBinding:
      position: 1

outputs:
  faprotax_output:
    label: Contains all FAPROTAX output
    type: Directory
    outputBinding:
      glob: $(inputs.output_directory_path)
  log_out:
    type: stdout
stdout: faprotax-tool.log

baseCommand:
- python
arguments:
- prefix: --out_collapsed
  position: 7
  valueFrom: faprotax_output/funct_otu_table.tsv
- prefix: -c
  position: 6
  valueFrom: '#'
- prefix: --out_sub_tables_dir
  position: 8
  valueFrom: faprotax_output
- prefix: -v
  position: 9
- prefix: --force
  position: 10
$schemas:
- http://edamontology.org/EDAM_1.16.owl
- https://schema.org/version/latest/schemaorg-current-http.rdf
s:copyrightHolder:
- name: Haris Zafeiropoulos
- url: https://hariszaf.github.io/
s:license: https://www.gnu.org/licenses/gpl-3.0.en.html
