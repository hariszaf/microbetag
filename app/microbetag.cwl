#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool


baseCommand: python3

inputs:
   parameters:
      type: File
      inputBinding:
         position: 1

outputs:
   example_out:
      type: stdout

stdout: output.txt






# baseCommand: python3

# inputs:
#    parameters:
#       type: File
#       inputBinding:
#          position: 1

#    otu_table: 
#       type: File
#       inputBinding:
#          position: 1
   
#    metadata:
#       type: File
#       inputBinding:
#          position: 1
   
#    edgelist:
#       type: File
#       inputBinding:
#          position: 1

# outputs:
#    example_out:
#       type: stdout

# stdout: output.txt

   


