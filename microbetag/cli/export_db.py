"""
Export the pre-calculations of a genome catalogue (the output folder of the microbetag-nf precalc workflow)
as the tables of a microbetag database (microbetag-export-db).

The output folder has one tab-separated file per table, without header, in the column order of
server/services/db/schema/custom_init.sql, and the seed sets the server computes seed complementarities from:

    genome_taxonomy.tsv           genomeTaxonomy:           genome, taxonomy
    phen_predictions.tsv          phenDB:                   genome, then a presence (YES/NO/N/A) and a score per trait
    unique_complements.tsv        uniqueComplements:        complement id, KEGG module, complement KOs, KOs of the path
    pathway_complementarity.tsv   pathwayComplementarity:   beneficiary, donor, complement ids (comma-separated)
    genome_ids.tsv                patricId2genomeId:        genome, genome (the catalogue's own identifiers)
    seed_scores.tsv               seedScores:               genome A, genome B, competition, cooperation
    seed_complements.tsv          SeedComplements:          beneficiary, donor, compounds (comma-separated)
    non_seed_sets.tsv             nonSeedSets:              genome, compounds (comma-separated)
    all_conf.json.gz, all_nonseeds.json.gz                   seed set confidences and non-seed sets per genome
"""

import argparse
import gzip
import json
import os
import shutil
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


def pathway_tables(pcompls_json: str) -> tuple[list[list], list[list]]:
    """
    From microbetag's pathway complementarities ({beneficiary: {donor: [[module, complement KOs, path KOs, url],
    ...]}}): the unique complements, each with an id, and the complement ids of each genome pair.
    """
    with open(pcompls_json) as f:
        pcompls = json.load(f)

    unique: dict[tuple[str, str, str], int] = {}
    pairs = []
    for beneficiary, donors in pcompls.items():
        for donor, complements in donors.items():
            if beneficiary == donor or not complements:
                continue
            ids = []
            for complement in complements:
                module, complement_kos, path_kos = complement[0], complement[1], complement[2]
                key = (
                    str(module).replace("md:", ""),
                    ";".join(complement_kos) if isinstance(complement_kos, list) else str(complement_kos),
                    ";".join(path_kos) if isinstance(path_kos, list) else str(path_kos),
                )
                if key not in unique:
                    unique[key] = len(unique) + 1
                ids.append(unique[key])
            pairs.append([beneficiary, donor, ",".join(map(str, dict.fromkeys(ids)))])

    unique_rows = [[cid, *key] for key, cid in unique.items()]
    return unique_rows, pairs


def seed_tables(seed_dir: str) -> tuple[list[list], list[list], list[list]]:
    """Seed scores, seed complements and non-seed sets from the seed_compl folder of precalc."""
    scores = pd.read_csv(os.path.join(seed_dir, "seed_scores.tsv"), sep="\t", dtype=str)
    if not scores.columns[0].startswith("node"):
        scores = pd.read_csv(os.path.join(seed_dir, "seed_scores.tsv"), sep="\t", dtype=str, header=None)
    score_rows = scores.iloc[:, :4].values.tolist()

    with open(os.path.join(seed_dir, "all_compls.json")) as f:
        compls = json.load(f)
    compl_rows = [
        [beneficiary, donor, ",".join(cpds)]
        for beneficiary, donors in compls.items()
        for donor, cpds in donors.items()
        if beneficiary != donor
    ]

    with open(os.path.join(seed_dir, "nonseeds.json")) as f:
        nonseeds = json.load(f)
    nonseed_rows = [[genome, ",".join(cpds)] for genome, cpds in nonseeds.items()]

    return score_rows, compl_rows, nonseed_rows


def gzip_copy(src: str, dest: str) -> None:
    with open(src, "rb") as f_in, gzip.open(dest, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def export(precalc: str, genome_taxonomy: str, outdir: str) -> dict[str, int]:
    """Write the database files; returns the number of rows of each table."""
    from microbetag.taxonomy import load_genome_taxonomy

    os.makedirs(outdir, exist_ok=True)
    counts = {}
    out = lambda name: os.path.join(outdir, name)

    taxonomy = load_genome_taxonomy(genome_taxonomy)
    counts["genomeTaxonomy"] = write_tsv(sorted(taxonomy.items()), out("genome_taxonomy.tsv"))

    predictions = os.path.join(precalc, "phenotrex", "predictions")
    if os.path.isdir(predictions):
        counts["phenDB"] = write_tsv(phen_table(predictions), out("phen_predictions.tsv"))

    pcompls = os.path.join(precalc, "pathway_compl", "pcompls.json")
    if os.path.exists(pcompls):
        unique, pairs = pathway_tables(pcompls)
        counts["uniqueComplements"]      = write_tsv(unique, out("unique_complements.tsv"))
        counts["pathwayComplementarity"] = write_tsv(pairs, out("pathway_complementarity.tsv"))

    seed_dir = os.path.join(precalc, "seed_compl")
    if os.path.isdir(seed_dir):
        scores, compls, nonseeds = seed_tables(seed_dir)
        counts["seedScores"]      = write_tsv(scores, out("seed_scores.tsv"))
        counts["SeedComplements"] = write_tsv(compls, out("seed_complements.tsv"))
        counts["nonSeedSets"]     = write_tsv(nonseeds, out("non_seed_sets.tsv"))
        gzip_copy(os.path.join(seed_dir, "confidence.json"), out("all_conf.json.gz"))
        gzip_copy(os.path.join(seed_dir, "nonseeds.json"), out("all_nonseeds.json.gz"))

    # The catalogue's genomes are their own "PATRIC" ids (the GTDB database maps GTDB genomes to PATRIC GEMs)
    genomes = set(taxonomy)
    counts["patricId2genomeId"] = write_tsv(sorted((g, g) for g in genomes), out("genome_ids.tsv"))

    missing = {k for k in ("phenDB", "pathwayComplementarity", "seedScores") if k not in counts}
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
    args = parser.parse_args(argv)

    for table, n in export(args.precalc, args.genome_taxonomy, args.outdir).items():
        print(f"{table:24s} {n} rows")


if __name__ == "__main__":
    main()
