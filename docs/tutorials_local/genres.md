---
title: Building GENREs
layout: default
parent: Additional tutorials
nav_order: 3
description: "how to build GENREs on microbetag"
---



# Building GENREs on `microbetag`


`microbetag` supports 2 ways to reconstruct GEMs based on the user's genomes/bins: 

1. using the [`modelseedpy`](https://github.com/ModelSEED/ModelSEEDpy) Python library
2. using the [`CarveMe`](https://carveme.readthedocs.io/en/latest/) tool 


In the first case, `modelseedpy` requires <a href="https://rast.nmpdr.org" target="_blank">RAST</a>-annotated genomes.
`microbetag` can do that on its own starting from your genome sequences; 
alternative, you may provide these to be used for the GEM reconstruction directly if you already have them
(either from previous `microbetag` runs or from other software).


```{note}
`modelseedpy` needs to establish a connection to the RAST server (`RastClient()`)
In some cases, based on the status of the RAST server, we have observed that time errors may occur. 
In this case, `microbetag` will exit and force a restart of its running on its own! 
Yet, it is a good practice to also check its status when the `modelseed` reconstruction step is running.
```


In the following paragraphs, we highlight how to go for different scenarios of GEMs reconstruction using different file types as initial starting points. 
One need to combine 2 parameters of the `config.yml` file to specify those scenarios: the `sc_input_type` where one specifies the file type and the `sequence_files_for_reconstructions` that points to the directory where the files to be used are located.



### using `modelseedpy` and your bins 

in this case, you have set 

- `sc_input_type` as `bins_fasta`, and 
- `sequence_files_for_reconstructions` is blank
- `genre_reconstruction_with` as `modelseedpy`

Then, `microbetag` will use [`RASTtk` programs](https://www.bv-brc.org/docs///cli_tutorial/rasttk_getting_started.html) to RAST annotate the original genomes/bins. 
In the `output_directory`, a folder called `reconstructions` has been built and in this case, 3 files for each genome/bin are now available:

- `.gto` and `.gto_2`: these are genome typed object, i.e. JSON files that are compatible with KBase. The `.gto_2` is a second genome typed object with all the RAST annotation data.
- `.faa` includes the same information as the `.gto_2` file, but we export the protein translations in `.fasta` format

```{note}
For our 7 genomes/bins this step may take about 1 hour depending on your computing system
```



### using `modelseedpy` and your already RAST annotated genomes

Assuming you already have the `.faa` files coming from the `rast-tk` package, you may use them directly by setting 

- `sc_input_type` as `proteins_faa`, and 
- `sequence_files_for_reconstructions` as the path to the folder with your `.faa` files
- `genre_reconstruction_with` as `modelseedpy`

In this case, `microbetag` will have to establish connections with the RAST client like before. 

```{note}
If your annotated genomes include the DNA sequences instead of the protein ones (`.fna` files) you may use them by setting the 
`sc_input_type` as `coding_regions`.
```


### using `carveme`

- `sc_input_type` as `bins_fasta`
- `sequence_files_for_reconstructions` is blank
- `genre_reconstruction_with` as `carveme`

In this case, under the `reconstructions` file, we have a `.tsv` file for each genome/bin with the findings of the `diamond` against the internal database of `carveme` with the BiGG reactions. 

|              |                 |       |     |      |     |     |      |    |      |            |        |
|:------------:|----------------:|------:|----:|-----:|----:|----:|-----:|---:|-----:|-----------:|-------:|
| bin_151.peg.3 |  iLJ478.TM0057 |  57.9 | 309 | 125  |  3  |  6  | 310  | 2  | 309  | 2.72e-128  |  369   |
| bin_151.peg.3 |  iLJ478.TM1063 |  55.9 | 311 | 130  |  3  |  7  | 310  | 3  | 313  | 5.36e-124  |  358   |

For a thorough description of each column, you may check this [here](https://github.com/bbuchfink/diamond_docs/blob/master/1%20Tutorial.MD).





