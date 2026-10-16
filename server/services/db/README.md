## Your own database, from your own genomes

Run the `precalc` workflow of the pipeline (`pipeline/`) on your genome catalogue, then export its output as the
files of a database, with the taxonomy of each genome (e.g. GTDB-Tk's `gtdbtk.*.summary.tsv`, or a tab-separated
genome and lineage):

```bash
microbetag-export-db --precalc <precalc output folder> --genome-taxonomy genome_taxonomy.tsv --outdir db_files
```

[`schema/microbetagDB.sql`](./schema/microbetagDB.sql) creates the tables and loads these files. MySQL runs it at its
first start, from `/docker-entrypoint-initdb.d/`, with the files in `/var/lib/mysql-files/`:

```bash
docker run -d --name microbetag-db \
    -e MYSQL_ROOT_PASSWORD=<password> \
    -v $PWD/schema/microbetagDB.sql:/docker-entrypoint-initdb.d/microbetagDB.sql:ro \
    -v $PWD/db_files:/var/lib/mysql-files \
    -p 3306:3306 mysql:8.4
```

The `/var/lib/mysql-files` mount needs to be writable: MySQL changes its owner at start (on Fedora and other
SELinux systems, add `:Z` to both mounts). Instead of matching taxa to genomes through NCBI taxonomy ids, as the
GTDB database does, a database of your own matches them by taxonomy (`genomeTaxonomy` table), to the genomes that
share their deepest rank.


## The GTDB database

## The database container

The idea is this 

- Dockerfile

- docker-compose.yml

    Gives credential for the database, as they are provided in the .env file

- .env

    The credentials to login to the database

- init.sql

    Describes the schema

- data/

    Symbolic link where we keep data mapped to the /var/lib/mysql-files


By running

    cd services/db/
    docker compose up

you will get a running container just with the database.


To have it as part of the 3-containers setup, see the main [README]()


## Luna specific

Since we cannot use for storage an account under `/home`, we keep data to be added in the database, under the `/var/lib/mysql-files`, and that is the actual directory we need to mount to our containers. 

However, for simplicity, one can make a symbolic link, pointing from that directory to the one where they are actually working, for example:

    sudo -i 
    cd /var/lib/mysql-files
    ln -s /home/u12324/microbetagApp/services/db/data .

> 🔴 To edit the `/var/lib/mysql-files` you need `sudo` rights. 


## Pathway complementarities 

We have a `.tar.gz` per genome with the actual complements:


    vsc35363@tier2-p-login-2 /staging/leuven/stg_00106/haris/microbetag/COMPLEMENTARITIES/to_compr_0 $ head complementarities_215221.tsv
    GCA_003207855.1	GCA_013282575.1	215221	2698676	md:M00121	K01845;K02492	K01885;K02492;K01845;K01698;K01749;K01719;K01599;K00228;K08973;K01772
    GCA_003207855.1	GCA_013282575.1	215221	2698676	md:M00846	K01845;K02492	K02492;K01845;K01698;K01749;K01719;K02302
    GCA_003207855.1	GCA_013282575.1	215221	2698676	md:M00868	K00231	K00643;K01698;K01749;K01719;K01599;K00228;K00231;K01772

Then, we have merged chunks of genomes in files with the prefix part_

    haris@pop-os:~$ head part_0_mapped_only_gtdb_genomes.tsv
    GCA_000222705.2	GCA_000005825.2	8409
    GCA_000222705.2	GCA_000005825.2	59264
    GCA_000222705.2	GCA_000005825.2	69574


where the third column corresponds to the complement id we have assigned to each unique complementarity found 



## Populate db

One can test/develop the database related parts on the host mysql, if available, but in the end of the day, we always need to make sure that we populate the database on the container

To do so, after firing the general `docker compose up` that inits all the 3 containers of the resource, we can run 

```bash
docker exec -it microbetagapp-db-1 bash
```

to jump into the container, and then:

```bash 
cd /docker-entrypoint-initdb.d/scripts/
./load_db.py
```

`load_db.py` connects with the credentials of the container (`MYSQL_USER`, `MYSQL_PASSWORD`, `MYSQL_DATABASE` of
`services/db/.env`), or with `USER_NAME`, `PASSWORD`, `HOST`, `PORT` and `DB_NAME` when it runs elsewhere.
`./load_db.py <env>` takes them from that environment of a `.env_dev.json` file instead.
It loads the GTDB database of our server; a database of your own genomes is loaded by `microbetagDB.sql`.



        "_comment": "USER_NAME in Luna used to be root"