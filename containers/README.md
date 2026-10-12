# Container images

The main image, with the microbetag package and the tools it calls, is built from the
[`Dockerfile`](../Dockerfile) at the root of the repository. It is used both stand-alone and by the
microbetag-nf Nextflow pipeline.

The images here hold single tools that need their own environment; the Nextflow pipeline runs them as
separate processes. Build them all from the **root** of the repository:

| Image | Dockerfile | Contents |
|---|---|---|
| `hariszaf/microbetag` | [`../Dockerfile`](../Dockerfile) | microbetag, Prodigal, HMMER, DIAMOND, FragGeneScan, CarveMe, manta, FlashWeave, phenotrex |
| `hariszaf/phenotrex` | `Dockerfile.phenotrex` | phenotrex 0.6 (same as `envs/mtg-phenotrex.yml`) |
| `hariszaf/modelseed` | `Dockerfile.modelseed` | ModelSEEDpy (same as `envs/mtg-modelseed.yml`) |
| `hariszaf/carveme` | `Dockerfile.carveme` | CarveMe and DIAMOND |
| `hariszaf/gapseq` | `Dockerfile.gapseq` | gapseq |

```bash
docker build -t hariszaf/microbetag:<version> .
docker build -f containers/Dockerfile.phenotrex -t hariszaf/phenotrex:0.6.0 .
```
