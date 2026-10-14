# microbetag
_annotating networks with phenotypic traits and potential metabolic interactions_

## About

`microbetag` annotates microbial co-occurrence networks:
  - on the **node** level: with phenotypic traits of the representing taxa
  - on the **edge** level: with potential metabolic interactions among the co-occurrying or mutual exclusive taxa.


The simplest way to use `microbetag` is through its Cytoscape App, called `MGG`.
You can [download `MGG`](https://apps.cytoscape.org/apps/mgg) directly from the Cytoscape ApppStore.

`MGG` allows both to run `microbetag` on-the-fly and visualize the returned annotated network.
However, using the on-the-fly version, means `microbetag` will try to match your taxa to their closest GTDB representative genomes.
Then, it will annotate your network based on pre-calculations of those GTDB genomes available on `microbetagDB`.

To annotate your network using your own genomes, run microbetag locally.
The steps that run external tools on the genomes (gene prediction, KEGG annotation, trait prediction,
genome-scale model reconstruction, FAPROTAX, network inference and clustering) are computationally
expensive; they run in the [**microbetag-nf**](pipeline) Nextflow pipeline,
on your computer or an HPC. The `microbetag` package then computes the pathway and seed complementarities
from their outputs and builds the annotated network, e.g. `mtag_net_2025-05-08_17-47.cx2`.
Load this file in Cytoscape to go through its annotations as in the on-the-fly version.

## Repository layout

The repository has three parts, one folder each:

```
microbetag/
├── microbetag/   the Python package: src/microbetag/, its tests, pyproject.toml and Dockerfile
├── pipeline/     the Nextflow pipeline (microbetag-nf)
├── server/       the on-the-fly server and its database (microbetagDB)
└── docs/         documentation (readthedocs)
```

| Part | What it does | Used by | Install / run |
|---|---|---|---|
| [**`microbetag/`**](microbetag) — the package | Computes pathway and seed complementarities from KEGG annotations and GEMs, and builds the annotated network (`.cx2`). It does not run external tools. Also exports pre-calculations as a database (`microbetag-export-db`) and matches taxonomies to genomes (`microbetag.taxonomy`). | you, the pipeline, the server | `pip install ./microbetag`, or the `hariszaf/microbetag` image; see [microbetag/README.md](microbetag/README.md) |
| [**`pipeline/`**](pipeline) — microbetag-nf | Runs the external tools in containers: Prodigal, kofam, phenotrex, CarveMe/gapseq/ModelSEEDpy, FAPROTAX, FlashWeave, manta, then the package. Workflows: `precalc` (your genomes → pre-calculations), `annotate` (→ an annotated network), `prep` (inputs for the on-the-fly version). | you, on a computer or an HPC | Nextflow and Docker/Singularity; see [pipeline/README.md](pipeline/README.md) |
| [**`server/`**](server) — on-the-fly server | The Flask app behind Cytoscape's `MGG`, and the MySQL microbetagDB it queries: GTDB genomes on the public server, or your own genomes after `precalc` and `microbetag-export-db`. | `MGG` | Docker Compose and systemd; see [server/README.md](server/README.md) |

How they fit together: the **pipeline** runs the tools and calls the **package**. The **package** builds the
annotated networks. The **server** serves the package's on-the-fly version from a database built from the
pipeline's `precalc` output.

If you already have the outputs of those tools, for example GEMs or a KEGG annotation of your genomes,
you can run `microbetag` on them directly.


For more on how to use `microbetag`, `MGG` and their features, feel free to have a look 
at [`microbetag`'s documentation page](https://microbetag.readthedocs.io/en/v1.0.3/).


## Installation: Running `microbetag`..

### Option 1: On-the-fly 

To use `microbetag` on-the-fly, all you have to do is to make sure you have Cytoscape on your computing system and add the `MGG` app.

* **You can [install Cytoscape](https://cytoscape.org/) from their website.**

* **After firing Cytoscape you can [add `MGG` by downloading it from the Cytoscape AppStore](https://apps.cytoscape.org/apps/mgg).**


### Option 2: The microbetag-nf pipeline

To go from your genomes to an annotated network, use the [microbetag-nf](pipeline)
pipeline. It runs the external tools in containers and then `microbetag`.


### Option 3: The microbetag package

`microbetag` needs Python 3.10 or later. From a clone of this repository:

```bash
git clone https://github.com/hariszaf/microbetag.git
cd microbetag/microbetag
pip install .
```

Add `".[otf]"` for the on-the-fly version, which queries the microbetag database, and `".[test]"` to run the tests.

To check your installation, run the tests from the same folder:

```bash
pip install ".[test]"
pytest
```

[`microbetag/tests/README.md`](./microbetag/tests/README.md) describes how to run single tests, or single microbetag steps on the test data.

#### In a container

The `hariszaf/microbetag` image holds the package and its commands:

```bash
docker pull hariszaf/microbetag:<version>
```

or build it from the `microbetag/` folder of this repository:

```bash
docker build -t hariszaf/microbetag:<version> .
```

Mount the folder with your configuration YAML and the files it points to at `/data`.
Paths in the YAML are relative to the YAML's folder, so a configuration that works locally also works in the container:

```bash
docker run --rm -v $PWD/my_study:/data hariszaf/microbetag:<version> microbetag --config /data/config.yml
```

For Singularity / Apptainer, convert the image once (no `sudo` needed) and run it with the same mount:

```bash
singularity pull microbetag.sif docker://hariszaf/microbetag:<version>
singularity exec -B $PWD/my_study:/data microbetag.sif microbetag --config /data/config.yml
```


## Further `microbetag` links

Overall, `microbetag` is a software ecosystem composed of a collection of modular code components:

- `microbetag-nf`: the Nextflow pipeline running the external tools, and preparing amplicon data for the on-the-fly version; in [`pipeline/`](pipeline)
- `microbetagDB` : the database of the server, in [`server/services/db`](server/services/db)
- `microbetag_prep`: replaced by the prep workflow of `microbetag-nf`; the former tool is on the `preprocess` branch of the [msysbio/microbetag GitHub repo](https://github.com/msysbio/microbetag/tree/preprocess)
- `microbetagApp`: individual [GitHub repo](https://github.com/msysbio/microbetagApp-public)
- `MGG` source code: individual [GitHub repo](https://github.com/hariszaf/MGG)
- pre-calculations: Seed and non-seed sets for GTDB representative genomes available in `microbetag`, their corresponding GEMs, models for `phenotrex`-based predictions and KEGG annotations of each of those genomes are available on [Zenodo](https://zenodo.org/records/15102937).
  
For potential issues when running `microbetag`, or ideas for further features, 
please feel free to join us on [`Matrix`](https://matrix.to/#/#microbetagcommunity:matrix.org).



## How to run `microbetag`

### On-the-fly

In case you are about to run `microbetag` using the on-the-fly option, you need to only provide you actual input data files, i.e.:
- an abundance file (mandatory)
- a co-occurrence network (optional)
- metadata file for inferring the co-occurrence network (optional)

For thorough descriptions on the format of each of these files, 
[continue on the `MGG` tutorials](https://microbetag.readthedocs.io/en/latest/basic_usage/mgg_totorials.html).


### Locally

`microbetag` needs a **configuration YAML file**. Each `microbetag` version comes with its own template,
which this command writes out:

```bash
microbetag --write-template config.yml
```

It points to:
- the **co-occurrence network**, as an edge list, and an **abundance table** or a sequence id to taxonomy map;
- for each annotation you want, the **output of the tool** it comes from:

| Annotation | Input | Produced by |
|---|---|---|
| Pathway complementarity | `ko_merged_file`: KEGG ORTHOLOGY terms per genome | kofam (hmmsearch) on the predicted genes |
| Seed complementarity | `genome_scale_models`: one GEM (.xml) per genome | CarveMe, gapseq or ModelSEEDpy |
| Phenotypic traits | `phenotrex_predictions` | phenotrex |
| FAPROTAX traits | `faprotax_sub_tables` | FAPROTAX |
| Network clusters | `prev_clustered_network` | manta |

microbetag-nf produces all of them. Switch off (`value: false`) the steps you have no input for.
Then run:

```bash
microbetag --config config.yml
```

[`microbetag/tests/data/microbetag`](./microbetag/tests/data/microbetag/) is a complete example, with its [`config.yml`](./microbetag/tests/data/microbetag/config.yml).

The steps are also available as separate commands (`microbetag-pcompl`, `microbetag-seed-sets`,
`microbetag-annotate`, ...), which the microbetag-nf pipeline uses; each has `--help`.
`microbetag` can also be used as a Python library, e.g. to get the seed complementarities of a set of GEMs
without a network (see [`microbetag/tests/unit/test_seed_compl.py`](./microbetag/tests/unit/test_seed_compl.py)).

## How microbetag gets seed sets

The seed complementarities of `microbetag` start from the **seed set** of each genome-scale metabolic model (GEM):
the compounds the organism has to take from its environment, because nothing in its metabolic network produces
them (Borenstein et al., [PNAS 2008](https://doi.org/10.1073/pnas.0712149105)). The rest of the compounds of the
network are its **non-seeds**: what it makes. microbetag follows the approach of
[PhyloMint](https://github.com/mgtools/PhyloMint) (Levy & Borenstein, [PNAS 2013](https://doi.org/10.1073/pnas.1300926110);
Lam et al., [PLoS Comput Biol 2020](https://doi.org/10.1371/journal.pcbi.1007951)), in
[`PhyloMint/lib/BuildGraphNetX.py`](./microbetag/src/microbetag/PhyloMint/lib/BuildGraphNetX.py).

Seeds are topological, so they do not depend on the medium the model was saved with: a seed is what the network
cannot make from anything, needed or not.

### 1. The compound graph of a GEM (`buildDG`)

The GEM (SBML, read with COBRApy) becomes a directed graph of compounds: an edge from each substrate to each product
of a reaction, in the directions the reaction's bounds allow (both, for a reversible one). Then:

- **one node per compound**: the forms of a compound in the compartments of the model (`ac_c`, `ac_e`, `ac_p`) are
  merged into one node (`ac`), so that seeds and non-seeds are per compound;
- **transports are left out**: a reaction with the same compounds on both sides, in other compartments (uniport,
  symport, antiport, ATP-driven), converts nothing, and would otherwise link co-transported compounds to each
  other (Ca²⁺ "made" from Na⁺);
- **every compound of the model is a node**, also one that is only transported or only used by the biomass (e.g.
  metal ions): one that nothing produces is then a source, as it should be;
- **the biomass reaction is a use**: its substrates get an edge to a sink node (`__biomass__`), which is neither a seed
  nor a non-seed;
- **currency metabolites** are left out, since they would link everything to everything: water, H⁺, ATP/ADP/AMP,
  phosphate and pyrophosphate, CO₂, O₂, NAD(P)(H), CoA, NH₄⁺, SO₄²⁻, Fe²⁺/Fe³⁺, ubiquinone and FAD
  (`CURRENCY` in `BuildGraphNetX.py`). The nutrients among them (NH₄⁺, SO₄²⁻, iron) get their status in a later
  round (section 2);
- **compounds made from currency metabolites alone** (ADP-ribose and NMN from NAD, cAMP and Ap4A from ATP,
  bicarbonate from CO₂) get an edge from a node standing for the currency pool (`__currency__`): every cell has
  these cofactors, so such compounds are non-seeds, not sources;
- exchange, demand and sink reactions, and reactions with zero bounds, are left out. Optionally (`check_blocked`),
  so are reactions that cannot carry flux on any medium (FVA with every exchange open); it is off by default, as it
  does not improve the seed sets of the models below.

### 2. Seeds and non-seeds (`getSeedSet`)

The graph is split into strongly connected components (SCCs: compounds that can be made from each other). An SCC
is a **source** if no edge enters it from outside. Then:

- a **seed** is a compound of a source SCC that also feeds the rest of the network (`require_downstream_use`, on by
  default: at least one edge out of the SCC, the biomass included). Its **confidence** is 1 / the size of its SCC;
- a source SCC that feeds nothing (a compound the model could import but never uses) is a **dead end**: neither a
  seed nor a non-seed, since nothing produces it;
- SCCs larger than `maxComponentSize` (2 in microbetag) are hubs, not sources;
- the **non-seeds** are the other compounds: what the network makes;
- **nutrient round**: NH₄⁺, SO₄²⁻ and iron (`NUTRIENT_CURRENCY`) are non-seeds if the network makes them from its
  seeds alone, e.g. NH₄⁺ from the deamination of an amino acid or amine, SO₄²⁻ from a sulfonate. This is network
  expansion over the model's reactions: a reaction fires once its substrates other than currency metabolites are
  available, starting from the seeds only, so that nitrogen, sulfur or iron cannot come from the cofactors (NH₄⁺
  from AMP). Otherwise they are seeds, with confidence 1. Iron released from an imported siderophore complex
  (ferric enterobactin, ferric dicitrate) does not count as made: the complex comes from outside, so iron stays a
  seed.

The seed and non-seed sets of each model are then translated to ModelSEED compound ids (BiGG models, through
MetaNetX).

### 3. Scores and complements of a pair (`CalculateIndexes.py`)

For a beneficiary A and a donor B:

- **competition** = the confidence-weighted share of A's seeds that are also B's seeds: both need them;
- **cooperation** = |A's seeds ∩ B's non-seeds| / |A's seeds ∩ (B's seeds ∪ B's non-seeds)|: of A's seeds that are in
  B's network, the share B makes;
- **seed complements** of A from B = A's seeds ∩ B's non-seeds: what B could provide A; by default, kept to the
  compounds of KEGG modules (the scores use the whole sets).

### Checking the seed sets

[`microbetag/scripts/check_seed_sets.py`](./microbetag/scripts/check_seed_sets.py) checks the seed sets of any GEM
against the model itself, with flux balance analysis:

- **disjoint**: no compound is both a seed and a non-seed;
- **completeness**: the *required* nutrients (the model does not grow without them even on the complete medium, so
  nothing in the network replaces them) are seeds;
- **precision**: seeds the model can make (with their own uptake closed, a demand reaction for them carries flux),
  which should not be seeds;
- biomass precursors reached by network expansion from the seeds.

```
cd microbetag
python scripts/check_seed_sets.py tests/data/seed_compl/iML1515.xml.gz tests/data/seed_compl/input_files/*.xml
```

| Model | Seeds | Required nutrients not seeds | Seeds the model can make |
|---|---|---|---|
| iML1515 (*E. coli* K-12, curated) | 126 | 0 of 9 | 2 (Fe²⁺, Fe³⁺) |
| `bin_101` (CarveMe, tests) | 55 | 3 of 11 | 0 |
| `GPB:bin_000083` (CarveMe, tests) | 41 | 10 of 28 | 1 (Fe³⁺) |

The seeds the models can make are iron, released from imported siderophore complexes: seeds by the definition above,
which flux balance analysis counts as made. In iML1515, NH₄⁺ is a non-seed (from dopamine, by an amine oxidase) and so
is SO₄²⁻ (from H₂S); iron is a seed.

With the former seed sets (cytosol only, one node per compartment form), iML1515 had 8 of its 9 required nutrients
(its metal ions and chloride) missing from its seeds, and 42 of its 127 seeds could be made by the model.

**What the checks can and cannot tell.** Seeds are not required nutrients by definition: a seed is what the network
cannot make from anything. The FBA checks are indicative, and their limits are those of the checks, not of the
inference. A nutrient a model requires but whose graph has a reaction producing it (whose other substrates, or
stoichiometry, the model cannot satisfy) is a non-seed, and shows as a "required nutrient not a seed": e.g. benzoate,
nicotinamide mononucleotide or tryptophan in the CarveMe models of the tests. A compound that is not required can
still be a seed (an alternative carbon source), and one required on a given medium need not be (glucose on M9 for
*E. coli*, which its network makes from other sugars).

**Known limits.** Two BiGG compounds can share a ModelSEED id (L-arabinitol and L-lactaldehyde, `cpd00334`), which can
then be both a seed and a non-seed after the translation. The currency metabolites other than NH₄⁺, SO₄²⁻ and iron are
never seeds.

**Tests**: [`test_seed_inference.py`](./microbetag/tests/unit/test_seed_inference.py) checks that seeds and non-seeds
are disjoint and that iML1515's required nutrients are seeds;
[`test_get_seeds.py`](./microbetag/tests/unit/test_get_seeds.py) checks the seed rules on small graphs. To draw the
graphs of the latter as well (needs `matplotlib`), from the `microbetag/` folder:

```
pytest tests/unit/test_get_seeds.py --draw --graph-dir tests/data/test_graphs
```

<!-- Issue where different findings from carveme using cplex and gurobi is mentioned:
https://github.com/cdanielmachado/carveme/issues/141#issuecomment-912309490 -->



## RTD

`microbetag` uses ReadTheDocs and Sphinx for its documentation (see `docs/`).

> **For contributors**
>
> To test changes locally first, run:
>
> Make sure you install `requirements.txt` on a conda environment: 
> 
>```
> cd docs/
> pip install -r requirements.txt
>```
>
> and then, you may run:
> ```
> sphinx-build -b html -d _build/doctrees -D language=en . _build/html -v
> ```
> from within the `docs/` directory.
> 

A new RTD is being produced every time a new _microbetag_ `tag` is released automatically.



## Cite

If you use `microbetag`, please consider citing us:

Zafeiropoulos H, Michail Delopoulos EI, Erega A, Schneider A, Geirnaert A, Morris J, Faust K.
microbetag: simplifying microbial network interpretation through annotation, enrichment tests and metabolic complementarity analysis.
*Genome Biology* 26, 292 (2025). DOI: [10.1186/s13059-025-03769-2](https://doi.org/10.1186/s13059-025-03769-2)


## Funding

This project is funded by: 
- the [3D' omics](https://www.3domics.eu) Horizon project (101000309).
- an [EMBO Short-Term Fellowship](https://www.embo.org/funding/fellowships-grants-and-career-support/scientific-exchange-grants/)

