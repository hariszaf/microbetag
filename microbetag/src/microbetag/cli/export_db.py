"""
Export the pre-calculations of a genome catalogue (the output folder of the microbetag-nf precalc workflow)
as the tables of a microbetag database (microbetag-export-db).

The output folder has one tab-separated file per table, without header, in the column order of
server/services/db/schema/microbetagDB.sql, and the seed sets the server computes seed complementarities from:

    genome_taxonomy.tsv           genomeTaxonomy:           genome, taxonomy
    phen_predictions.tsv          phenDB:                   genome, then a presence (YES/NO/N/A) and a score per trait
    unique_complements.tsv        uniqueComplements:        id, KEGG module, missing KOs, KOs of the path: each
                                                            alternative way to complete a module, once
    genome_alternatives.tsv       genomeAlternatives:       genome, ids of the alternatives of its incomplete modules
    genome_module_kos.tsv         genomeModuleKos:          genome, KEGG module, KOs of the module it has
    seed_sets.tsv                 seedSets:                 GEM, its seeds with their confidence (cpd:1.0,cpd:0.5)
    non_seed_sets.tsv             nonSeedSets:              GEM, its non-seeds (comma-separated)
    genome_ids.tsv                patricId2genomeId:        genome, id of its GEM (by default the genome itself)
    genome_ncbi_tax_ids.tsv       genome2taxNcbiId:         NCBI Taxonomy id, genome (empty without --ncbi-tax-ids)

The server computes the pathway and seed complementarities of the pairs of a network from these, so nothing is
stored per pair. From precalc: phenotrex/predictions, pathway_compl/alternatives.json (or .json.gz),
hmmsearch/ko_merged.txt (or .txt.gz) and seed_compl (confidence.json, nonseeds.json, or .json.gz).
"""

import argparse
import gzip
import json
import os
import re
from collections.abc import Iterable

import pandas as pd

# phenotrex model (prediction file name) -> phenDB column, in the column order of phenDB
PHEN_TRAITS: list[tuple[str, str]] = [
    ("ac", "aceticAcid"), ("AOB", "aob"), ("a_saccharolytic", "aSaccharolytic"), ("auto_co2", "autoCo2"),
    ("butanol", "butanol"), ("butyric_acid", "butyricAcid"), ("d_glucose", "dGlucose"),
    ("d_lactic_acid", "dLacticAcid"), ("ethanol", "ethanol"), ("fermentative", "fermentative"),
    ("fixing_n2", "fixingN2"), ("formic_acid", "formicAcid"), ("halophilic", "halophilic"),
    ("hydrogen", "hydrogen"), ("indole", "indole"), ("isobut", "isobutyricAcid"),
    ("isovaleric_acid", "isovalericAcid"), ("l_lactic_acid", "lLacticAcid"), ("methanotroph", "methanotroph"),
    ("NOB", "NOB"), ("non_fermentative", "nonFermentative"), ("phototrophy", "phototrophy"),
    ("psychrophilic", "psychrophilic"), ("r_acetoin", "rAcetoin"), ("saccharolytic", "saccharolytic"),
    ("succinic_acid", "succinicAcid"), ("sulfate_reducer", "sulfateReducer"), ("symbiont", "symbiont"),
    ("T3SS", "T3SS"), ("T6SS", "T6SS"), ("thermophylic", "thermophilic"), ("anaerobe", "anaerobe"),
    ("aerobe", "aerobe"),
]

FASTA_SUFFIXES = (".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna", ".faa")


def genome_id(name: str) -> str:
    """Genome identifier of a genome file name, e.g. bin_101.fa -> bin_101."""
    name = os.path.basename(str(name).strip())
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def write_tsv(rows: Iterable[Iterable], path: str) -> int:
    """
    Write rows as a tab-separated file without header; returns the number of rows.
    Missing values (None) are written as \\N, which MySQL's LOAD DATA reads as NULL.
    """
    n = 0
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join("\\N" if v is None else str(v) for v in row) + "\n")
            n += 1
    return n


def phen_table(predictions_dir: str) -> list[list]:
    """
    One row per genome: the genome, then the presence and the score of each trait, in phenDB column order.
    Traits without a prediction for the genome get N/A and no score.
    """
    traits: dict[str, dict[str, tuple[str, str]]] = {}
    genomes = set()
    for model, _ in PHEN_TRAITS:
        candidates = [f"{model}.tsv", f"{model}.prediction.tsv", f"{model}.predictions"]
        path = next((os.path.join(predictions_dir, c) for c in candidates
                     if os.path.exists(os.path.join(predictions_dir, c))), None)
        if path is None:
            continue
        # keep_default_na=False: phenotrex's "N/A" is a value, not a missing one
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
        df.columns = ["identifier", "presence", "confidence"][: len(df.columns)]
        traits[model] = {
            genome_id(r.identifier): (r.presence, r.confidence) for r in df.itertuples(index=False)
        }
        genomes.update(traits[model])

    rows = []
    for genome in sorted(genomes):
        row = [genome]
        for model, _ in PHEN_TRAITS:
            presence, confidence = traits.get(model, {}).get(genome, ("N/A", None))
            row += [presence, confidence]
        rows.append(row)
    return rows


KO = re.compile(r"^K\d{5}$")


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _first_existing(*paths: str):
    return next((p for p in paths if os.path.exists(p)), None)


def genome_kos(ko_merged: str) -> dict[str, set[str]]:
    """The KOs of each genome, from precalc's ko_merged.txt (genome, then the KO and the gene, in any order)."""
    kos: dict[str, set[str]] = {}
    with _open(ko_merged) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            ko = next((x for x in fields[1:] if KO.match(x)), None)
            if ko:
                kos.setdefault(genome_id(fields[0]), set()).add(ko)
    return kos


def pathway_tables(alternatives_json: str, ko_merged: str, module_kos_file: str) -> tuple[list, list, list]:
    """
    From precalc's alternatives (genome -> incomplete module -> alternative path -> its missing KOs) and KOs: each
    alternative once, with an id, the alternatives of each genome, and the KOs each genome has per module.
    """
    with _open(alternatives_json) as f:
        alternatives = json.load(f)

    unique: dict[tuple[str, str, str], int] = {}
    genome_rows = []
    for genome, modules in alternatives.items():
        ids = []
        for module, alts in modules.items():
            for path, missing in alts.items():
                path_kos = path if isinstance(path, list) else re.findall(r"K\d{5}", path)
                key = (module.replace("md:", ""), ";".join(sorted(missing)), ";".join(path_kos))
                ids.append(unique.setdefault(key, len(unique) + 1))
        genome_rows.append([genome_id(genome), ",".join(map(str, sorted(set(ids))))])

    module_kos = pd.read_csv(module_kos_file, sep="\t", header=None, names=["module", "ko"])
    module_kos["ko"] = module_kos["ko"].str.replace("ko:", "", regex=False)
    by_module = module_kos.groupby("module")["ko"].apply(set).to_dict()
    ko_rows = []
    for genome, kos in sorted(genome_kos(ko_merged).items()):
        for module, mkos in sorted(by_module.items()):
            present = sorted(kos & mkos)
            if present:
                ko_rows.append([genome, module.replace("md:", ""), ";".join(present)])

    unique_rows = [[cid, *key] for key, cid in unique.items()]
    return unique_rows, genome_rows, ko_rows


def two_columns(path: str) -> list[list[str]]:
    """Rows of a two-column tab-separated file."""
    with _open(path) as f:
        return [line.rstrip("\n").split("\t")[:2] for line in f if line.strip()]


def export(precalc: str, genome_taxonomy: str, outdir: str, gem_ids: str = None,
           ncbi_tax_ids: str = None) -> dict[str, int]:
    """
    Write the database files; returns the number of rows of each table.

    gem_ids: tab-separated genome and the id of its GEM, when they differ (the GTDB database: PATRIC ids)
    ncbi_tax_ids: tab-separated NCBI Taxonomy id and genome, for taxa given as NCBI Taxonomy ids
    """
    from microbetag.helpers import MappingPaths
    from microbetag.taxonomy import load_genome_taxonomy

    os.makedirs(outdir, exist_ok=True)
    counts = {}
    out = lambda name: os.path.join(outdir, name)

    taxonomy = load_genome_taxonomy(genome_taxonomy)
    counts["genomeTaxonomy"] = write_tsv(sorted(taxonomy.items()), out("genome_taxonomy.tsv"))

    predictions = os.path.join(precalc, "phenotrex", "predictions")
    if os.path.isdir(predictions):
        counts["phenDB"] = write_tsv(phen_table(predictions), out("phen_predictions.tsv"))

    pc = os.path.join(precalc, "pathway_compl")
    alternatives = _first_existing(os.path.join(pc, "alternatives.json"), os.path.join(pc, "alternatives.json.gz"))
    hmm = os.path.join(precalc, "hmmsearch")
    ko_merged = _first_existing(os.path.join(hmm, "ko_merged.txt"), os.path.join(hmm, "ko_merged.txt.gz"))
    if alternatives and ko_merged:
        module_kos = os.path.join(MappingPaths().kegg_mappings, "kegg_terms_per_module.tsv")
        unique, genome_alts, genome_module_kos = pathway_tables(alternatives, ko_merged, module_kos)
        counts["uniqueComplements"]  = write_tsv(unique, out("unique_complements.tsv"))
        counts["genomeAlternatives"] = write_tsv(genome_alts, out("genome_alternatives.tsv"))
        counts["genomeModuleKos"]    = write_tsv(genome_module_kos, out("genome_module_kos.tsv"))

    sc = os.path.join(precalc, "seed_compl")
    confidence = _first_existing(os.path.join(sc, "confidence.json"), os.path.join(sc, "confidence.json.gz"))
    nonseeds = _first_existing(os.path.join(sc, "nonseeds.json"), os.path.join(sc, "nonseeds.json.gz"))
    if confidence and nonseeds:
        with _open(confidence) as f:
            conf = json.load(f)
        with _open(nonseeds) as f:
            nons = json.load(f)
        counts["seedSets"] = write_tsv(
            ([gem, ",".join(f"{c}:{v:g}" for c, v in seeds.items())] for gem, seeds in sorted(conf.items())),
            out("seed_sets.tsv"))
        counts["nonSeedSets"] = write_tsv(
            ([gem, ",".join(cpds)] for gem, cpds in sorted(nons.items())), out("non_seed_sets.tsv"))

    # The GEM of each genome: the genome's own id, unless given
    gems = two_columns(gem_ids) if gem_ids else sorted((g, g) for g in taxonomy)
    counts["patricId2genomeId"] = write_tsv(gems, out("genome_ids.tsv"))
    counts["genome2taxNcbiId"]  = write_tsv(two_columns(ncbi_tax_ids) if ncbi_tax_ids else [],
                                            out("genome_ncbi_tax_ids.tsv"))

    missing = {k for k in ("phenDB", "genomeAlternatives", "seedSets") if k not in counts}
    if missing:
        print(f"Not in {precalc}, so not exported: {', '.join(sorted(missing))}")
    return counts


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-export-db", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--precalc", required=True, help="Output folder of the precalc workflow")
    parser.add_argument("--genome-taxonomy", required=True,
                        help="Tab-separated genome and lineage of each genome of the catalogue, e.g. GTDB-Tk's summary")
    parser.add_argument("--outdir", required=True, help="Folder for the database files")
    parser.add_argument("--gem-ids", help="Tab-separated genome and the id of its GEM, when they differ")
    parser.add_argument("--ncbi-tax-ids", help="Tab-separated NCBI Taxonomy id and genome, if known")
    args = parser.parse_args(argv)

    counts = export(args.precalc, args.genome_taxonomy, args.outdir, args.gem_ids, args.ncbi_tax_ids)
    for table, n in counts.items():
        print(f"{table:24s} {n} rows")


if __name__ == "__main__":
    main()
