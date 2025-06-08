---
title: On your own genomes
layout: default
parent: Additional tutorials
nav_order: 3
description: "an example case of how to run microbetag using your own bins/MAGs"
---


# .. using custom genomes


This is the main versioning scheme for `microbetag`. 
The version of this tutorial is in line with the one of the documentation page.

<div style="display: flex; gap: 10px;">
   <a href="https://hub.docker.com/r/hariszaf/microbetag" class="btn-green" target="_blank"> Docker image </a>
   <a href="https://github.com/hariszaf/microbetag/tree/user-bins/tests/dev_io_microbetag" class="btn-purple" target="_blank"> Tutorial files </a>
</div>



```{note}
This tutorial is for users that have some basic experience working on a terminal. 
 
Contrary to previous cases, this scenario is not performed from within the CytoscapeApp.

The user needs to run `microbetag` first on their computing environment (personal computer, HPC etc.), 
and then load the returned annotated network to Cytoscape.
```

<!-- 
You may find the input files we are using for this tutorial in the `user-bins` branch of `microbetag`s GitHub repo, 
under the [`tests/dev_io_microbetag` folder](https://github.com/hariszaf/microbetag/tree/user-bins/tests/dev_io_microbetag).
-->

In the [Cytoscape App tutorial](../tutorials_otf/abd_only.md), 
our sequences were already taxonomically assigned before running `microbetag`,
and their taxonomies were mapped to representative GTDB genomes.
`microbetag` then used these genomes for the annotation steps.


However, in case of shotgun metagenomics one may end up with their own bins/ 
while further refinement of the latter can lead to Metagenome-Assembled Genomes (MAGs).
Also, one may already use genomes for the communities under study from earlier studies.
In this case, such local genomes can be used directly for the annotation steps of `microbetag`. 
Yet, this requires computing resources and time much higher than those that our web-server can support.

<!-- 
that you may get from our 
[GitHub repo](https://github.com/hariszaf/microbetag/blob/user-bins/tests/dev_io_microbetag/config.yml) 
or download it directly from [here][1] 
-->

In this tutorial, we will use a very short number of bins (7) to showcase the various steps `microbetag` implements. 
Yet, it still gets more than a couple of hours to go through all the different steps on a personal computer.
In our experience, memory (RAM) requirements should not be a challenge; 
memory would be an issue only with really large networks.
`microbetag` is more often than not thread-limited, i.e. 
it needs computing power to go through the annotation steps. 


```{important}
**INPUT FILES USED IN THIS TUTORIAL**

 A complete example of running `microbetag` locally using the `modelseedpy` library for GEM reconstruction with all the intermediate files produced can be found in the 
 [`dev_io_microbetag`](https://github.com/hariszaf/microbetag/tree/user-bins/tests/dev_io_microbetag) folder of the `user-bins` 
 branch on the GitHub repo.
 
In the initial run, there are only 3 input files:
 - the <a href="{{ site.url }}/microbetag/download/local/config.yml" download="config.yml"><code>config.yml</code></a> file; allows you to set all the relative parameters for `microbetag` to run
<!-- [`config.yml`][1] file; allows you to set all the relative parameters for `microbetag` to run -->
 - an **abundance table** (following the format of the Cytoscape app tutorial) called [`thirty_Samples.tsv`][2], and
 - its corresponding edge list ([`edgelist.csv`][3])

 **Remember!**

 The **config** file is **mandatory**. 
```

```{hint}
This tutorial primarily covers the case where your starting input is an **abundance table** 
along with a list of genomes or bins.
Alternatively, you may also run `microbetag` using again the list of genomes but this time, 
a **precomputed network** file and a **two-column taxonomy file**, 
where the first column contains sequence IDs and the second contains their corresponding taxonomies.

For such a case, you may check an example on microbetag's 
<a href="https://github.com/msysbio/microbetag/tree/develop/test_data/test_microbetag" target="_blank">GitHub</a>.
```




## Input and `config.yml` files

The `config.yml` file is rather important as it is the one that allows you to set your `microbetag` run.
A number of the parameters there correspond to tools that are invoked,
while others have to do with alternative routes that `microbetag` can follow for the annotation of the network.
Read **carefully** the `description` of each argument before setting a value. 
Here, we highlight some of them. 


- `abundance_table_file`: path to your abundance table; the abundance table needs to follow the instructions 
  for any abundance table to be used with `microbetag`, i.e., sequence identifier in the first column, 
  sample names in the first row and a 7-level taxonomy in the last column; of course, 
  you may provide the output of the[ `microbetag` preprocessing step](../tutorials_otf/prep.md) as an abundance table. 

- `sc_input_type`: This is a **key parameter** for running `microbetag` locally; 
  based on whether you already have annotated your genomes (either using other software or from previous runs of `microbetag`) 
  you can use different input files as the starting point for getting the seed complementarities. 
  The `sequence_files_for_reconstructions` parameter is strongly related to this. 
  
  For example, if you have already GEMs reconstructed based on your genomes, you may set this to `sc_input_type` to `models` and then, 
  provide the folder name with your GEMs in the `sequence_files_for_reconstructions` parameter (e.g. `my_xmls`). 
  Likewise, if you do not have GEMs, but you already have RAST annotations, you may set `sc_input_type` to `proteins_faa` 
  and give the path to those in the `sequence_files_for_reconstructions` parameter.

- `seed_complementarity`: since this is the most time and resource consuming step, the user may choose not to go for it. 
  By setting this to `Fasle`, none of the steps for GEMs reconstruction or seed complementarity inference will be performed.

- `flashweave_args`: all the arguments under this umbrella term are related to how `FlashWeave` will perform, 
  check on the [FAQs](../faq.md#when-to-enable-the-sensitive-and-heterogeneous-arguments) but also the 
  <a href="https://github.com/meringlab/FlashWeave.jl" target="_blank">FlashWeave GitHub repo</a> for more.



```{warning}
Make sure that the template configuration file you are using is compatible with 
the version of `microbetag` you have installed.
```

```{hint}
A YAML configuration file is required **only** in case you need to run the whole *microbetag* pipeline.

If you need to run indipendent tasks, you can either build pseudo :class:`microbetag.Config` classes
or even use directly the `microbetag` features you wish to. 

The example cases we conver under the <a href="https://github.com/msysbio/microbetag/tree/develop/tests" target="_blank">`tests` folder on GitHub</a> they come along with their input/output data you will find on the 
<a href="https://github.com/msysbio/microbetag/tree/develop/test_data" target="_blank">`test_data`</a> folder.
Have a look on the <a href="https://github.com/msysbio/microbetag/blob/develop/test_data/README.md" target="_blank">`README`</a>
for the cases covered and feel free to advise their corresponding pseudo-config classes, or their corresponding
YAML configuration files.
```



## Output 

In the `config.yml` file, you can specify the `output_directory`.
Here we discuss the folders and the files you will find under the `output_directory`.
We dot not always follow the order with which the files are generated. 

### The annotated `.cx2` network file

The main output file of `microbetag` can be found in the `output_directory` you set in your `config.yml` file; 
the `microbetag`-annotated network called `mtag_net_<timestamp>.cx2`.
This is the file you need to [load in your Cytoscape](../tutorials_otf/prep.md) and then, a
after enabling the MGG visual style and the MGG results panel, you can investigate your annotated network! 
This file is in 
[`.cx2`](https://cytoscape.org/cx/cx2/specification/cytoscape-exchange-format-specification-(version-2)/){target="_blank"} 
format. 



### FAPROTAX
A folder called `faprotax` is made where there is a subfolder, called `sub_tables` 
and a file with the sum of the abundances of the taxa found with a specific process in each sample, called `functional_otu_table.tsv`. 
In the `sub_tables` folder, a file for each process is available mentioning the genomes/bins found 
related with the process under study and their relative abundance per sample. 

For example, the `aerobic_nitrite_oxidation.txt` looks like: 

| record |  seqId |   sample1 | sample2 | sample4 | sample5 | sample6 | sample7 | sample8 | sample9 | sample10 | sample11 | sample12 | sample13 |
| :----:|:-------:|:---------:|:-------:|:--------:|:-------:|:------:|:--------:|:------:|:-------:|:--------:|:--------:|:---------:|:-------:|
| d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Nitrobacter;s__Nitrobacter sp001896955 | bin_32  | 43 |  10 | 56  | 73 | 9 | 58 | 54 | 46 | 9  | 40 | 42 | 81 |
| d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Nitrobacter;s__Nitrobacter sp001897285 | bin_223 | 47 |  87 | 69  | 64 | 25| 95 | 40 | 71 | 16 | 78 | 52 | 40 |



### genome-based trait predictions

- `train.genotype` file: this is the output of the `phenotrex` program annotating your genomes/bins with COG families using the latest 

- `predictions` folder: within this folder, a file with the predictions for presence/absence of each trait in each of our genomes along with a confidence score are available. For example, the `symbiont.prediction.tsv ` file, in our test case looks like: 

|# Trait: Symbiont|             |             |   
|:---------------:|:-----------:|:-----------:|
|Identifier	      |Trait present|	Confidence|
|bin_101.fa	      | NO	        | 0.643      |
|bin_151.fa	      | NO	        |  0.6395    |
|bin_19.fa	      | NO          |  0.869     |
|bin_38.fa        | NO	        |  0.8678    |
|bin_41.fa	      | NO	        |  0.7842    |
|bin_45.fa	      | YES	        |  0.7954   |
|bin_48.fa	      | NO	        |  0.8545   |



### ORFs

`microbetag` invokes `prodigal` to extract Open Reading Frames (ORFs). 
It creates a folder called `ORFs` in the `output_directory` and for each genome/bin it returns 3 files: 
- `.gbk`: Genbank-like format (for more information check [here](https://www.insdc.org/submitting-standards/feature-table/))
- `.faa`: the reading frames as aminoacid sequences
- `.ffn`: the reading frames as nucleic acid sequences


```{note}
**SKIP THE ORFs PREDICTION (`prodigal`) STEP**
 
If you have already calculated the ORFs of your genomes before start using `microbetag`, you can create a folder within your `output_directory` called `ORFs` and move there all your `.faa` files.
 This way, `microbetag` will be using those instead of running `prodigal`.

 You `.faa` files should look like this:


    >c_000000001749_1 # 2 # 913 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.479
    DDSKIHQLGWDAFQAGTKVAKEEGLYGAGQDLLSDAFSGNVKGLGPAVAELSFEERPSEP
    FLFFMADKTEPGAYNLPFYLSYADPMYNPGLMLSPKMGKGFVFTVMDVENTENDRIIELT
    TPEDIYDLACLLRDNGRFVVESIRSAKTGETTAVCSTTRLNKIAGEYVGKDDPVALARVQ

```


### KEGG annotations 

`microbetag` makes use of the `hmmsearch` tool and the `kofam_database` profiles to check which KOs are present in each of your genomes.
In the `output_directory`, `microbetag` creates a folder called `KEGG_annotations` and there it builds a folder called `hmmout`, where it keeps all the 24.728 `.hmmout` files for each genome.  
Once all the `.hmmout` files are there for all the genomes/bins under study, `microbetag` builds a file called `ko_merged.txt` based on the DiTing implementation, that looks like this:

| bin_id	|    contig_id                       |	ko_term     |
|:---------:|:----------------------------------:|:------------:|
| bin_41	| SCN18_26_2_15_R1_F_scaffold_115_57 |	K07586      |
| bin_48	| SCN18_26_2_15_R4_B_scaffold_93_80	 |  K08086      |
| bin_41	| SCN18_26_2_15_R1_F_scaffold_206_63 |	K03503      |

This file is the main component for `microbetag` to proceed with the pathway complementarity step. 


```{tip}
**SKIP THE KEGG ANNOTATION (`HMMSEARCH`) STEP**
 
If you have already hmm profiles either from analysis before using `microbetag` or from previous `microbetag` runs of your genomes, you can create a folder called `hmmout` within the `KEGG_annotations` folder of your `output_directory` and move all the `.hmmout` profiles of your bins there. 
An `.hmmout` file would looks like:

    root@8649bd465c24:/data/microbetag_local/KEGG_annotations/hmmout# more K00005.bin_101.hmmout 
    #                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
    ># target name        accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
    #------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
    c_000000006615_9     -          K00005               -            1.2e-98  327.2   0.0   1.4e-98  327.0   0.0   1.0   1   0   0   1   1   1   1 # 14663 # 15772 # 1 # ID=74_9;partial=00;start_type=TTG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;gc_cont=0.449
    #
    # Program:         hmmsearch
    # Version:         3.4 (Aug 2023)
    # Pipeline mode:   SEARCH
    # Query file:      /microbetag/microbetagDB/ref-dbs/kofam_database/profiles/K00005.hmm
    # Target file:     /data/microbetag_local/ORFs/bin_101.faa
    # Option settings: hmmsearch -o /dev/null --tblout /data/microbetag_local/KEGG_annotations/hmmout/K00005.bin_101.hmmout -T 324.27 --cpu 1 /microbetag/microbetagDB/ref-dbs/kofam_database/profiles/K00005.hmm /data/microbetag_local/ORFs/bin_101.faa 
    # Current dir:     /microbetag
    # Date:            Mon Aug  5 11:06:08 2024
    # [ok]

 If you already have the `ko_merged.txt` file, you can only add a copy of it in the `KEGG_annotations` folder (the `hmmout` files are not necessary in this case) and `microbetag` will use this directly skipping the `hmmsearch` step. 
```




```{hint}
**COMPUTING TIME, RESOURCES AND STORAGE**

Running `microbetag` locally using your own genomes/bins/MAGs can take significant computing time and resources.
In this tutorial, we use only a short number of bins that has almost no biological significance.
What we want to accomplish here is to make sure that you can run `microbetag` locally.
Indicatively, using a Linux machine and allocating 2 CPUs it took more than 1 hour for the KO annotation of just those 7 bins.
Running the complete workflow could get up to 6-7 hours based on the approach you will choose to reconstruct your GEMs. 
 
`carveme` is much more robust in running smoothly and faster since it does not require a RAST connection (see following paragraph). 
```




### GEMs already available

In this case, you may use your GEMs directly for the seed complementarities inference by setting:

- `sc_input_type` as `models`
- `sequence_files_for_reconstructions` pointing to directory with the `.xml` files
- `genre_reconstruction_with` can be left blank or any value; it will not be considered




### PhyloMInt post process

After `microbetag` performs `PhyloMInt`, it runs a step to post-process the seed and the non seed sets as initially returned by:
- removing compounds from seed sets that are related to environmental metabolites that can be produced in several ways within the cell.
- removing from non seed sets compounds that cannot be produced in any other way than from entering the cell from the environment.

The numbers of this post process are recorded in the `log.tsv` file that can be found under the `seed_complementarity` folder that looks like this:

|model_id | environmental_initial_seeds | non_environmental_initial_seeds | total_initial_seeds | updated_seeds | initial_non_seeds | updated_non_seeds |
|:-------:|:---------------------------:|:-------------------------------:|:-------------------:|:-------------:|:-----------------:|:-----------------:|
| bin_151 |            173              |                57               |            230      |      230      |        879        |       879         |
| bin_38  |            177              |                54               |            231      |      231      |       1211        |      1211         |
| bin_101 |            152              |                51               |            203      |      203      |        845        |       845         |

This post process step is necessary since we use a complete medium to gapfill the model. 

```{note}
One may come with alternative procedures on how to gap fill in terms of minimising the missing potential cross-feedings, but at the same time not over-predicting such cases.
```


[1]:../_static/download/local/config.yml
[2]:../_static/download/local/thirty_Samples.tsv
[3]:../_static//download/local/edgelist.csv

