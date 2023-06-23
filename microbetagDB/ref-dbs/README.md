# Data and databases used in `microbetag`

```bash
for NCBI in */ ; do cd $NCBI ; sed -i 's/ko://g' *.json  ; cd .. ; done
```


Currently, `microbetag` exploits genomes from the following databases:

* KEGG
* GTDB
* MGnify catalogues

| Resource | unique ncbi ids added |
| :------: | :-------------------- |
|   GTDB   | 14,818                |
|   KEGG   | 3,984                 |
|  MGnify  | 1,044                 |



Short chunk of code to concatenate all the annotated genomes under the `all_genome_modules` folder.

```
import os
import json

directory = os.path.join( os.getcwd(),  "all_genome_modules")

json_files = []

for root, dirs, files in os.walk(directory):
    for file in files:
        file_path = os.path.join(root, file)
        json_files.append(file_path)


data_dict = {}

for file_path in json_files:
    with open(file_path, "r") as file:
        ncbi_id = file_path.split("/")[-2]
        file_data = json.load(file)
        if ncbi_id not in data_dict:
            data_dict[ncbi_id] = {}
            data_dict[ncbi_id][file_path] = file_data
        else:
            data_dict[ncbi_id][file_path] = file_data

out_file = open("ko_per_mod_in_genome_per_ncbiId.json", "w")
json.dump(data_dict, out_file)
```

NCBI Taxonomy Ids to neglect:

- 297314 > uncultured Lachnospiraceae bacterium
- 1898203 > Lachnospiraceae bacterium
 

