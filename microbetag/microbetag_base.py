#!/usr/bin/env python3

# Information about the microbetag tool
version      = "v0.1.0"
license      = ("GPLv3",)
packages     = ["microbetag"]
description  = "a microbial co-occurrence network annotator"
author       = "Haris Zafeiropoulos"
author_email = "haris.zafeiropoulos@kuleuven.be"
name         = "microbetag"

"""
microbetag is a microbial network annotator that exploits several software and databases to 
annotate microbial, co-occurence networks. 
Functional traits like whether a species is methanotroph, fermentative etc are assigned to 
each node of the network. 

For the associations that include 2 taxa of the species/strain level, microbetag also performs a 
pathway complementarity step; species are assigned to their corresponding GTDB representative 
genomes and using those genomes, complementarity for KEGG modules are investigated. 

Finally, microbetag supports a series of metabolic models covering the GTDB representative genomes. 
Metabolic analysis methods such as FBA and flux sampling are performed to highlight potential metabolic 
interactions among the taxa. 
"""


import os
from microbetag.scripts.variables import *
from microbetag.scripts.logger import *
from microbetag.scripts.db_functions import *
from microbetag.scripts.utils import *


def main():

    """
    Assure the output directory
    """
    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)

    """
    Setting logging; the logger variable comes from the logger.py script
    """
    # Using FileHandler writing log to file
    logfile = os.path.join(OUT_DIR, "log.txt")
    fh      = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    # Using StreamHandler writing to console
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)

    # Using error
    eh = logging.StreamHandler()
    eh.setLevel(logging.ERROR)
    eh.setFormatter(formatter)

    # Add the two Handlers
    logger.addHandler(ch)
    logger.addHandler(fh)
    logger.addHandler(eh)

    """
    Welcome message and arguments values
    """
    logging.info("Hello microbe-fun! microbetag is about to start!")
    logging.info("Your command was: {}".format(" ".join(sys.argv)))

    """
    STEP 1: OTU table preprocess 
    """
    if OTU_TABLE: 

        otu_table = is_tab_separated(OTU_TABLE, TAX_COL)
        logging.info("STEP 1: Assign NCBI Tax Id and GTDB reference genomes".center(80, "*"))
        if not EDGE_LIST:

            logging.info("The user has not provided an edge list. microbetag will build one using FlashWeaeve.")
            if not os.path.exists(FLASHWEAVE_OUTPUT_DIR):
                os.mkdir(FLASHWEAVE_OUTPUT_DIR)

            ext = ensure_flashweave_format(otu_table, TAX_COL, OTU_COL)

        else:
            ext = otu_table.copy()
            ext["microbetag_id"] = otu_table[OTU_COL]

        # Map taxonomies to ontology ids
        logging.info("Get the NCBI Taxonomy id for those OTUs that have been assigned either at the species, the genus or the family level.")
        try:
            otu_taxid_level_repr_genome, repr_genomes_present = map_otu_to_ncbi_tax_level_and_id(ext, TAX_COL, OTU_COL)
        except:
            logging.error("microbetag was not able to map your table's taxonomies to NCBI and GTDB ids. Check on how you describe your table. Also check whether you have set the edge_list parameter in the config file; if there is not an edge list, leave it blank.")
            sys.exit(0)

        otu_taxid_level_repr_genome.to_csv( os.path.join( FLASHWEAVE_OUTPUT_DIR, "ncbi_parsed_otu.csv" ), "\t")

    """
    STEP 2: Get co-occurrence network
    """
    if not EDGE_LIST:
        """
        Run FlashWeave
        """
        logging.info("STEP 2: Get co-occurrence network".center(80, "*"))

        test_julia = os.system( "julia -v" )
        if test_julia != 0:
            logging.info("Julia is not present in the OS. Please intall or run microbetag as a Docker image.")
            sys.exit(0)

        flashweave_params = [
            "julia", FLASHWEAVE_SCRIPT, FLASHWEAVE_OUTPUT_DIR, FLASHWEAVE_TMP_INPUT
        ]
        flashweave_command = " ".join(flashweave_params)
        if os.system(flashweave_command) != 0:
            logging.error("FlashWeave was not able to perform. Check whether installed. You can either add FlashWeave or run microbetag as a Docker image. Also, check if issue with your OTU/ASV column.")
            sys.exit(0)

        # Taxa pairs as NCBI Tax ids
        logging.info("Map your edge list to NCBI Tax ids and keep only associations that both correspond to a such.")
        edge_list = edge_list_of_ncbi_ids(FLASHWEAVE_EDGELIST, otu_taxid_level_repr_genome)


    """
    STEP 3: PhenDB
    """
    logging.info("STEP 3: PhenDB ".center(80, "*"))
    if PHEN_DB: 

        # Check connectivity to db
        if not check_connection_to_db(USER_NAME, PASSWORD, HOST, DB_NAME):
            logging.error("microbetag cannot link to the microbetag db. Check host server, user name and password credentials.")
            sys.exit(0)

        # Get phen traits for each GTDB genome present on your table
        feats = get_column_names()
        feats.insert(1, "NCBI_ID")
        feats.insert(1, "Species")
        traits = []
        for gtdb_id in set(repr_genomes_present):
            check = True
            try:
                r = get_phenDB_traits(gtdb_id)
            except:
                check = False
            if not check:
                if gtdb_id[:3] == "GCA":
                    try:
                        r = get_phenDB_traits(gtdb_id.replace("GCA", "GCF"))
                        check = True
                    except:
                        pass
                else:
                    try:
                        r = get_phenDB_traits(gtdb_id.replace("GCF", "GCA"))
                    except:
                        pass
            if check:
                tr = otu_taxid_level_repr_genome.loc[otu_taxid_level_repr_genome["gtdb_gen_repr"] == gtdb_id, ["Species", "ncbi_tax_id"]]
                sp = tr.iloc[0,0]
                ncbiId = str(int(tr.iloc[0,1]))
                r.insert(1, ncbiId)
                r.insert(1, sp)
                traits.append(r)
            else:
                print("Genome :", gtdb_id, "is not present in the phenDB version of microbetag.")

        # Save phen traits as a .tsv file
        if not os.path.exists(PHEN_OUTPUT_DIR):
            os.mkdir(PHEN_OUTPUT_DIR)
        outfile = os.path.join(PHEN_OUTPUT_DIR, "phen_traits.tsv")
        export_phen_traits_to_file(column_names = feats, rows = traits, filename=outfile)

    """
    STEP 4: FAPROTAX
    """
    logging.info("STEP 4: FAPROTAX database oriented analaysis".center(80, "*"))
    if OTU_TABLE: 

        if not os.path.exists(FAPROTAX_OUTPUT_DIR):
            os.mkdir(FAPROTAX_OUTPUT_DIR)

        faprotax_params = [
            "python3", FAPROTAX_SCRIPT,
            "-i",      OTU_TABLE,
            "-o",      FAPROTAX_FUNCT_TABLE,
            "-g",      FAPROTAX_DB,
            "-c", '"' + COM_CHAR + '"',
            "-d", '"' +  TAX_COL + '"',
            "-v",
            "--force",
            "-s",      FAPROTAX_SUB_TABLES,
        ]

        if COM_HEAD:
            faprotax_params = faprotax_params + ["--column_names_are_in", COM_HEAD]

        faprotax_command = " ".join(faprotax_params)

        # Run FAPROTAX 
        # In the sub tables files, in column 1 we have the taxonomy and in column 2 the OTU ID. 
        if os.system(faprotax_command) != 0:
            logging.error("\nSomething went wrong when running the BugBase analysis! Check how you describe your input table. Also, please make sure you have set the column_names_are_in parameter properly.")
            sys.exit(0)

        path_to_subtables = os.path.join(BASE, FAPROTAX_SUB_TABLES)
        otu_faprotax_assignments = otu_faprotax_functions_assignment(path_to_subtables)
        print(otu_faprotax_assignments)
        """
        [TODO] Link OTU ID or the taxonomy to the GTDB ID if such!
        """

    """
    STEP 5: BugBase
    """
    logging.info("STEP 5: BugBase database oriented analaysis".center(80, "*"))
    if BUGBASE and OTU_TABLE: 

        # Make a copy of the otu table without the taxonomy column 
        f = open(OTU_TABLE, "r")
        g = open(BUGBASE_TMP, "w")
        for line in f:
            g.write("\t".join(line.split("\t")[:-1]) + "\n")

        bugbase_params = [
            "Rscript", BUGBASE_SCRIPT, 
            "-i", BUGBASE_TMP,
            "-o", BUGBASE_OUTPUT, 
            "-a", 
        ]

        if METADATA_FILE:
            bugbase_params = bugbase_params + ["-m", METADATA_FILE]

        bugbase_command = " ".join(bugbase_params)
        print(bugbase_command) 
        print(">>>", BUGBASE_PATH)

        # Run BugBase
        logging.info(["Command to run: ", bugbase_command])
        if os.system(bugbase_command) != 0:
            logging.error("\nSomething went wrong when running the BugBase analysis!")
            sys.exit(0)

        """
        [TODO]: - Parse the bugbase/otu_contributions/contributing_otus.txt to assign features in the OTUs 
        """
        os.remove(BUGBASE_TMP)

        # for tr_file in os.listdir(BUGBASE_OUTPUT)


    sys.exit(0)

#     """
#     STEP 6: PATHWAY COMPLEMENTARITY
#     """
#     logging.info("STEP 6: Pathway complementarity".center(80, "*"))
#     if PATHWAY_COMPLEMENTARITY: 

#         """
#         Remember to change the path from kegg_genomes to all_genomes once the latter is ready 
#         """
#         set_of_ncbi_ids_with_available_genomes = set(os.listdir("ref-dbs/kegg_genomes/"))

#         if EDGE_LIST:

#             edge_list = EDGE_LIST.copy()

#         edgelist_to_pc = edge_list.to_dict(orient="records")
#         """
#         Example:
#         [{'node_a': 'microbetag_17', 'ncbi_tax_id_node_a': 77133.0, 'gtdb_gen_repr_node_a': 'GCA_903925685.1', 'ncbi_tax_level_node_a': 'mspecies', 'node_b': 'microbetag_21', 'ncbi_tax_id_node_b': 136703.0, 'gtdb_gen_repr_node_b': nan, 'ncbi_tax_level_node_b': 'species'}, 
#         {'node_a': 'microbetag_17', 'ncbi_tax_id_node_a': 77133.0, 'gtdb_gen_repr_node_a': 'GCA_903925685.1', 'ncbi_tax_level_node_a': 'mspecies', 'node_b': 'microbetag_74', 'ncbi_tax_id_node_b': 77133.0, 'gtdb_gen_repr_node_b': 'GCA_903925685.1', 'ncbi_tax_level_node_b': 'mspecies'}... 
#         """

#         for pair in edgelist_to_pc: 

#             if pair["ncbi_tax_level_node_a"] == pair["ncbi_tax_level_node_b"] == "mspecies":

#                 taxon_a = str(pair["ncbi_tax_id_node_a"]).split(".")[0]
#                 taxon_b = str(pair["ncbi_tax_id_node_b"]).split(".")[0]

#                 print(pair["gtdb_gen_repr_node_a"], pair["gtdb_gen_repr_node_b"])

#                 if taxon_a in set_of_ncbi_ids_with_available_genomes and taxon_b in set_of_ncbi_ids_with_available_genomes: 
#                     print("hello friend")
#                     print(pathway_complementarity(taxon_a, taxon_b))

#     return True


#     """
#     STEP 7: Net Cooperate
#     """
#     logging.info("STEP 7: NetCooperate between species/strains paired nodes".center(80, "*"))
#     # if NETCOOPERATE: 



main()

