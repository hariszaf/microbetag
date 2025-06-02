


docker run --rm -it  \
--volume=/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/datasets/dev_test:/data \
--volume=/home/luna.kuleuven.be/u0156635/github_repos/KU/microbetag/ext_data/kofam_database/:/microbetag/microbetagDB/ref-dbs/kofam_database/ \
--volume=/home/luna.kuleuven.be/u0156635/github_repos/KU/microbetag/gurobi.lic:/opt/gurobi/gurobi.lic:ro \
--entrypoint /bin/bash  microbetag:v1.0.7 



# FOR MGG
# ==========
# script to fix on  MGG for the panel to have all annotations (those "_" are not there)
# check whether the same issue affects the enrcichment analysis

# src/main/java/be/kuleuven/mgG/internal/utils/Mutils.java


docker run --rm -it  --volume=/home/luna.kuleuven.be/u0156635/Documents/projects/microbetag/datasets/validation/:/data --volume=/home/luna.kuleuven.be/u0156635/github_repos/KU/microbetag/ext_data/kofam_database/:/microbetag/microbetagDB/ref-dbs/kofam_database/ --volume=/home/luna.kuleuven.be/u0156635/github_repos/KU/microbetag/gurobi.lic:/opt/gurobi/gurobi.lic:ro --entrypoint /bin/bash  microbetag:v1.0.7 


I should have seen a 2nd info. Percentage of mapped metabolites: 98.68421052631578
