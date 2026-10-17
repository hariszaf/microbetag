"""
The GTDB database's inputs in the layout of the pipeline's precalc output, so that it is exported with
microbetag-export-db as any catalogue:

    precalc/phenotrex/predictions/<trait>.tsv   from gtdb_phen_predictions_oxy.tsv (phenDB's columns)
    precalc/pathway_compl/alternatives.jsonl.gz from the KO files (all_genome_modules/), one genome per line
    precalc/hmmsearch/ko_merged.txt.gz          the same KOs, as genome, KO, gene ("-")
    precalc/seed_compl/{confidence,nonseeds}.json.gz   the seed sets of the GEMs (PATRIC ids), from run_seed_sets.py
    genome_taxonomy.tsv, gem_ids.tsv, ncbi_tax_ids.tsv

Genomes: those of gtdb2patricIds.tsv (GTDB r207 representatives with a PATRIC GEM). Inputs, in /l: the files of
the former GTDB database (gtdb2patricIds.tsv, update_genomeId2ncbiTaxId.tsv, gtdb_phen_predictions_oxy.tsv); in /e:
the KO files of the genomes (all_genome_modules/, Zenodo record 8316356), the seed sets, and GTDB r207's
bac120/ar53 metadata. Then:

    python build_gtdb_precalc.py <workers>
    microbetag-export-db --precalc precalc --genome-taxonomy genome_taxonomy.tsv --outdir db_files \
        --gem-ids gem_ids.tsv --ncbi-tax-ids ncbi_tax_ids.tsv
"""
import gzip, json, os, sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from microbetag.cli.export_db import PHEN_TRAITS
from microbetag.helpers import MappingPaths
from microbetag.pathway_complementarity import process_bin

L, E = "/l", "/e"
OUT = f"{E}/gtdb_db/precalc"
num = lambda g: g.split("_", 1)[1].split(".")[0]
m = MappingPaths()
ref = pd.read_csv(os.path.join(m.kegg_mappings, "kegg_terms_per_module.tsv"), sep="\t", header=None, names=["module", "ko"])
MODULES = sorted(ref["module"].unique())
MO_MAP = json.load(open(os.path.join(m.kegg_mappings, "module_definition_map.json")))
STRUCTURALS = ["md:M00144", "md:M00149", "md:M00151", "md:M00152", "md:M00154", "md:M00155", "md:M00153",
               "md:M00156", "md:M00158", "md:M00160"]


def one(args):
    genome, path = args
    present = {mod: sorted(set(v.values())) for mod, v in json.load(open(path)).items()}
    _, alts = process_bin(genome, {mod: present.get(mod, []) for mod in MODULES}, MO_MAP, STRUCTURALS)
    kos = sorted({k for v in present.values() for k in v})
    return genome, {mod: {p: sorted(g) for p, g in a.items()} for mod, a in alts.items()}, kos


def main(workers):
    for d in ("phenotrex/predictions", "pathway_compl", "hmmsearch", "seed_compl"):
        os.makedirs(f"{OUT}/{d}", exist_ok=True)
    g2p = [l.rstrip("\n").split("\t") for l in open(f"{L}/gtdb2patricIds.tsv") if l.strip()]
    genomes = {g for g, _ in g2p}
    by_num = {num(g): g for g in genomes}

    # Genomes, their GEMs and NCBI Taxonomy ids
    with open(f"{E}/gtdb_db/gem_ids.tsv", "w") as f:
        f.writelines(f"{g}\t{p.strip()}\n" for g, p in g2p)
    with open(f"{E}/gtdb_db/ncbi_tax_ids.tsv", "w") as f:
        for l in open(f"{L}/update_genomeId2ncbiTaxId.tsv"):
            t, g = l.rstrip("\n").split("\t")[:2]
            if g.startswith(("GCA_", "GCF_")):
                f.write(f"{t}\t{g}\n")
    lineage = {}
    for meta in ("bac120_metadata_r207.tsv", "ar53_metadata_r207.tsv"):
        df = pd.read_csv(f"{E}/../Documents/unclassified/{meta}", sep="\t", usecols=["accession", "gtdb_taxonomy"],
                         dtype=str)
        for acc, tax in df.itertuples(index=False):
            lineage[num(acc[3:])] = tax
    with open(f"{E}/gtdb_db/genome_taxonomy.tsv", "w") as f:
        f.writelines(f"{g}\t{lineage[num(g)]}\n" for g in sorted(genomes) if num(g) in lineage)
    print(f"genomes {len(genomes)}, with a GTDB r207 lineage {sum(num(g) in lineage for g in genomes)}", flush=True)

    # phenotrex predictions, one file per trait
    phen = pd.read_csv(f"{L}/gtdb_phen_predictions_oxy.tsv", sep="\t", header=None, dtype=str, keep_default_na=False)
    for i, (model, _) in enumerate(PHEN_TRAITS):
        with open(f"{OUT}/phenotrex/predictions/{model}.tsv", "w") as f:
            f.write(f"# Trait: {model}\nIdentifier\tTrait present\tConfidence\n")
            for g, presence, conf in zip(phen[0], phen[1 + 2 * i], phen[2 + 2 * i]):
                f.write(f"{g}\t{presence}\t{conf}\n")
    print(f"phenotypes of {len(phen)} genomes", flush=True)

    # Seed sets of the GEMs (run_seed_sets.py: one JSON per GEM), merged
    confidence, nonseeds = {}, {}
    sets = f"{E}/gtdb_seed_sets/sets"
    for name in sorted(os.listdir(sets)):
        if name.endswith(".json"):
            d = json.load(open(os.path.join(sets, name)))
            confidence[d["genome"]], nonseeds[d["genome"]] = d["confidence"], d["nonseeds"]
    for data, dst in ((confidence, "confidence"), (nonseeds, "nonseeds")):
        with gzip.open(f"{OUT}/seed_compl/{dst}.json.gz", "wt") as f:
            json.dump(data, f)
    print(f"seed sets of {len(confidence)} GEMs", flush=True)

    # KOs and alternatives, streamed
    todo = []
    for d, _, fs in os.walk(f"{E}/gtdb_db/all_genome_modules"):
        for name in fs:
            n = num(name.replace("_kos_related_to_mos.json", ""))
            if n in by_num:
                todo.append((by_num[n], os.path.join(d, name)))
    todo.sort()
    with gzip.open(f"{OUT}/pathway_compl/alternatives.jsonl.gz", "wt") as fa, \
            gzip.open(f"{OUT}/hmmsearch/ko_merged.txt.gz", "wt") as fk, ProcessPoolExecutor(workers) as ex:
        fk.write("bin_id\tko_term\tgene\n")
        for i, (genome, alts, kos) in enumerate(ex.map(one, todo, chunksize=20), 1):
            fa.write(json.dumps({genome: alts}) + "\n")
            fk.writelines(f"{genome}\t{k}\t-\n" for k in kos)
            if i % 2000 == 0:
                print(f"{i} / {len(todo)}", flush=True)
    print(f"KOs and alternatives of {len(todo)} genomes", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]))
