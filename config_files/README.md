# Configuration files

Configuration templates per microbetag version. Copy the one of your version and set the paths of your files.

| Version | Template |
|---|---|
| 2.0 | [`config_v200.yml`](./config_v200.yml): microbetag works on the outputs of the external tools, which run in [microbetag-nf](https://github.com/msysbio/microbetag-nf) |
| 1.0.4 | [`config_v104.yml`](./config_v104.yml): microbetag also runs the external tools |

In 2.0, each step that used to run a tool takes that tool's output instead:

| Step | 1.0.4 | 2.0 |
|---|---|---|
| Network | inferred with FlashWeave if `edge_list` was missing | `edge_list` is required |
| Pathway complementarity | kofam annotation of `bins_fasta` with `kofam_database` | `ko_merged_file` |
| Seed complementarity | GEMs built from `bins_fasta` (`sc_input_type`, `genre_reconstruction_with`, ...) | `genome_scale_models` |
| Phenotypic traits | phenotrex on `bins_fasta` | `phenotrex_predictions` |
| FAPROTAX | run on the abundance table | `faprotax_sub_tables` |
| Network clustering | manta | `prev_clustered_network` |

1.0.4 also used `pc_percentage`, `prev_calc_path_compl` and `precalulations_only`, which microbetag
did not read; the keys it reads are `pc_percent`, `pc_file` and `precalculations_only`.
