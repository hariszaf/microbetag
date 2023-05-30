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
