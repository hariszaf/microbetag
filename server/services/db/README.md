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
./load_db.py <env>
```
`<env>` should be an environment mentioned in the `.env_dev.json` file. 



        "_comment": "USER_NAME in Luna used to be root"