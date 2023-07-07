# NCBI Tax Id - Genomes available - KO terms found links

The `microbetag_ko_db.json` file is the *core* mapping file for the pathway complementarity step. 
It links the NCBI Taxonomy Ids supported by `microbetag` with their corresponding genomes that have been retrieved.
For example, here is an entry of this `.json` file:
```
    "1655": {
        "tax-level": "species",
        "resources": [
            "kegg_genomes",
            "mgnify_catalogues"
        ],
        "links": {
            "kegg_genomes": {
                "ane": {
                    "inner-link": "ref-dbs/kegg_genomes/1655/ane_kos_related_to_mos.json",
                    "outer-link": "https://rest.kegg.jp/link/ko/ane"
                }
            },
            "mgnify_catalogues": {
                "MGYG000299101": {
                    "inner-link": "ref-dbs/mgnify_catalogues/1655/MGYG000299101_kos_related_to_mos.json",
                    "outer-link": "http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-oral/v1.0/species_catalogue/MGYG0002991/MGYG000299101/genome"
                }
            }
        }
    }
```





