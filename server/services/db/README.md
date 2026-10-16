# The microbetag database

Every microbetag database, that of our server (GTDB) and one of your own genomes, has the same tables
([`schema/microbetagDB.sql`](./schema/microbetagDB.sql)), and is built the same way: the pre-calculations of a
genome catalogue (the pipeline's `precalc` workflow), exported as one file per table with `microbetag-export-db`,
then loaded by MySQL. The server reads, for each network, the rows of its genomes, and computes their pathway and
seed complementarities; nothing is stored per pair of genomes.

| Table | One row per | Content |
|---|---|---|
| `genomeTaxonomy` | genome | its lineage; taxa are matched to the genomes sharing their deepest rank (`CATALOGUE=custom`) |
| `genome2taxNcbiId` | genome | its NCBI Taxonomy id, for taxa given as NCBI ids (the GTDB database); empty otherwise |
| `patricId2genomeId` | genome | the id of its GEM (PATRIC ids for the GTDB database; the genome itself otherwise) |
| `phenDB` | genome | phenotrex predictions: presence and score of each trait |
| `uniqueComplements` | alternative | a way to complete a KEGG module: the module, its missing KOs, the KOs of its path |
| `genomeAlternatives` | genome | the alternatives of its incomplete modules |
| `genomeModuleKos` | genome and module | the KOs of the module it has (for when it is the donor) |
| `seedSets` | GEM | its seeds, with their confidence |
| `nonSeedSets` | GEM | its non-seeds |


## Your own database

Run the `precalc` workflow of the pipeline (`pipeline/`) on your genome catalogue, then export its output, with
the taxonomy of each genome (e.g. GTDB-Tk's `gtdbtk.*.summary.tsv`, or a tab-separated genome and lineage):

```bash
microbetag-export-db --precalc <precalc output folder> --genome-taxonomy genome_taxonomy.tsv --outdir db_files
```

MySQL runs `schema/microbetagDB.sql` at its first start, from `/docker-entrypoint-initdb.d/`, with the files in
`/var/lib/mysql-files/`:

```bash
docker run -d --name microbetag-db \
    -e MYSQL_ROOT_PASSWORD=<password> -e MYSQL_DATABASE=microbetagDB \
    -e MYSQL_USER=<user> -e MYSQL_PASSWORD=<password> \
    -v microbetag-db-data:/var/lib/mysql \
    -v $PWD/schema/microbetagDB.sql:/docker-entrypoint-initdb.d/microbetagDB.sql:ro \
    -v $PWD/db_files:/var/lib/mysql-files \
    mysql:8.4
```

or, with the server, `docker compose up` from `server/` (`services/db/data` is the folder of the files, see
[docker-compose.yml](../../docker-compose.yml)). The `/var/lib/mysql-files` mount needs to be writable, as MySQL
changes its owner at start, and on a Linux file system (not NTFS or FAT); on Fedora and other SELinux systems,
add `:Z` to the mounts.

MySQL loads the files only when its data folder (`/var/lib/mysql`, the `microbetag-db-data` / `mysql-data` volume)
is empty. To load new files, remove that volume (`docker compose down -v`, or `docker volume rm`) and start
again.


## The GTDB database of our server

Built the same way, from the GTDB r207 representative genomes with a ModelSEED GEM (33,756 GEMs; KEGG annotations
for 16,098 of them), with the scripts of [`gtdb/`](./gtdb):

1. `run_seed_sets.py`: the seed sets of the GEMs (Zenodo record 23368237);
2. `build_gtdb_precalc.py`: the inputs in the layout of `precalc`'s output, from the seed sets, the KEGG
   annotations (Zenodo record 8316356), the phenotrex predictions and the genomes' NCBI and PATRIC ids of the
   former database, and GTDB r207's metadata;
3. `microbetag-export-db` with `--gem-ids` (the PATRIC GEM of each genome) and `--ncbi-tax-ids`.

The exported files are about 550 MB; loaded, about 0.5 GB on disk, and MySQL uses about 0.6 GB of memory.
