import os
import re
import sys
import ast
import json
import math
import time
import cobra
import shutil
import random
import itertools
import pandas as pd
import multiprocessing
from modelseedpy import MSBuilder, MSGenome


class SetEncoder(json.JSONEncoder):
    """
    Custom JSON encoder that handles serialization of Python sets.

    This encoder extends the functionality of the standard JSONEncoder to support
    serializing Python sets. JSON does not have a native representation for sets,
    so this encoder converts sets to lists before serializing them.

    Usage:
        When serializing data to JSON using json.dump() or json.dumps(), specify
        cls=SetEncoder to use this custom encoder.

    References:
        - json.JSONEncoder: https://docs.python.org/3/library/json.html#json.JSONEncoder
    """
    def default(self, obj):
        """
        Override the default method of JSONEncoder to handle serialization of sets.
        Notes:
            If the object is a set, it is converted to a list before serialization.
            Otherwise, the default behavior of JSONEncoder.default() is used.        Notes:
            If the object is a set, it is converted to a list before serialization.
            Otherwise, the default behavior of JSONEncoder.default() is used.
        """
        if isinstance(obj, set):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


def get_files_with_suffixes(directory, suffixes):
    """

    """
    matching_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if any(file.endswith(suffix) for suffix in suffixes):
                matching_files.append(os.path.join(root, file))
    return matching_files


def run_prodigal(fasta, basename, outdir):
    """
    Function to predict ORFs using Prodigal.

    fasta:
    basename:
    outdir:
    """
    faa_file = os.path.join(outdir, basename + '.faa')
    cmd_para = [
                'prodigal', '-q',
                '-i', fasta,
                '-p', 'meta',
                '-a', faa_file,
                '-d', os.path.join(outdir, basename + '.ffn'),
                '-o', os.path.join(outdir, basename + '.gbk')
                ]
    cmd = ' '.join(cmd_para)
    if os.path.exists(faa_file):
        print("ORFs already predicted for bin:", basename)
    else:
        print("ORFs to be predicted for bin:", basename)
        try:
            os.system(cmd)
        except:
            print("Something wrong with prodigal annotation!")


def kegg_annotation(faa, basename, out_dir, db_dir, ko_dic, threads):
    """
    Function to perform KEGG annotation.
    The function invokes hmmsearch.

    faa:
    basename:
    out_dir:
    db_dir:
    ko_dic:
    threads:
    """
    print('KEGG annotation for {}'.format(basename))
    paras = []
    for knum, info in ko_dic.items():

        output = os.path.join(out_dir, knum + '.' + str(basename) + '.hmmout')
        hmm_db = os.path.join(db_dir, 'profiles', knum + '.hmm')

        if os.path.exists(output):
            continue

        if not os.path.exists(hmm_db):
            continue

        if info[1] == 'full':
            threshold_method = '-T'
            outtype = '--tblout'

        elif info[1] == 'domain':
            threshold_method = '--domT'
            outtype = '--domtblout'

        elif info[1] == 'custom':
            threshold_method = '-E'
            outtype = '--tblout'

        paras.append((threshold_method, info[0], outtype, output, hmm_db, faa))

    print("Number of KEGG processes to be performed:", str(len(paras)))

    process = multiprocessing.Pool(threads)
    process.map(hmmsearch, paras)


def hmmsearch(paras):
    """
    Function to invoke hmmsearch software.
    """
    (threshold_method, threshold, outtype, output, hmm_db, faa) = paras
    cmd_para = [
        'hmmsearch',
        threshold_method, threshold,
        '--cpu', '1',
        '-o /dev/null',
        outtype, output,
        hmm_db,
        faa
    ]
    cmd = ' '.join(cmd_para)
    try:
        os.system(cmd)
    except:
        print("Something wrong with KEGG hmmsearch!")


def ko_list_parser(ko_list):
    """
    parse ko_list file into a dict object - based on DiTing

    :param ko_list: path of the file ko_list
    :return: a dictionary mapping knum to threshold and score_type
    :rtype: dict
    """
    ko_dic = {}  # { knum : [threshold, score_type] }
    with open(ko_list) as fi:
        next(fi)  # skip the first line (header)
        for line in fi:
            knum, threshold, score_type = line.split('\t')[0:3]
            if threshold == '-':
                continue
            else:
                ko_dic[knum] = [threshold, score_type]
    return ko_dic


def merge_ko(hmmout_dir, output):
    """
    Parses the KO<>.<bin>.hmmout files produced by the kegg_annotation() function
    to create a single 3-column file (output) with the bin_id, the corresponding conting and the KO that wa mapped to it.
    The function then returns a dictionary with the bin ids as the keys and the set of KOs found to each as the value.

    hmmout_dir: path to the .hmmout files
    output: path/filename to save the output file
    """
    if not os.path.exists(output):
        with open(output, 'w') as fo:
            fo.write('bin_id\tcontig_id\tko_term\n')
        for hmmout_file in os.listdir(hmmout_dir):
            if hmmout_file.endswith('.hmmout'):
                kobasename = hmmout_file.rsplit('.', 1)[0]
                basename = kobasename.split('.', 1)[1]
                hmmout_file_path = os.path.join(hmmout_dir, hmmout_file)
                with open(hmmout_file_path, 'r') as fi:
                    for line in fi:
                        if not line.startswith('#'):
                            gene_id, _ = line.split()[0:2]  # under _ the accession
                            lines = line.split()
                            if re.match(r'[0-9]+$', lines[2]):
                                k_number = lines[3]
                            else:
                                k_number = lines[2]
                            with open(output, 'a') as fo:
                                fo.write(basename + '\t' + gene_id + '\t' + k_number + '\n')

    df = pd.read_csv(output, sep="\t")

    bins_kos = df.groupby('bin_id')['ko_term'].apply(set).to_dict()

    # Pivot the DataFrame to have 'kegg_id' as rows and 'bin_id' as columns
    unique_combinations = df.drop_duplicates().copy()
    unique_combinations.loc[:, 'presence'] = 1
    pivot_df = unique_combinations.pivot_table(index='ko_term', columns='bin_id', values='presence', fill_value=0)

    return bins_kos, pivot_df  # keep one


def flatten(list_of_lists):
   """
   This function takes a list of lists and flattens it until it returns a list with
   all the components of the initial one.
   """
   if len(list_of_lists) == 0:
      return list_of_lists
   if isinstance(list_of_lists[0], list):
      return flatten(list_of_lists[0]) + flatten(list_of_lists[1:])
   return list_of_lists[:1] + flatten(list_of_lists[1:])


def build_kegg_url(kegg_map, clean_path, missing_kos):
   """
   Build url to colorify the related to the module kegg map based on the KO terms
   of the beneficiary (pink) and those it gets from the donor (green)
   """
   # Load the dictionary with the kegg modules and their corresponding maps
   color_mapp_base_url = "https://www.kegg.jp/kegg-bin/show_pathway?"
   present_kos_color   = "%09%23EAD1DC/"
   complemet_kos_color = "%09%2300A898/"

   # Make a url pointing at a colored kegg map based on what's on the beneficiary's genome and what it gets as complement from the donor
   beneficiarys_kos = ""
   complements_kos  = ""
   for ko_term in clean_path:
      if ko_term not in missing_kos:
         beneficiarys_kos = "".join([beneficiarys_kos, ko_term, present_kos_color])
      else:
         complements_kos = "".join([beneficiarys_kos, ko_term, complemet_kos_color])
   try:
      """ we do the try as in some rare cases, the module might not have a related map"""
      url_ko_map_colored = "".join([color_mapp_base_url, kegg_map,  "/", beneficiarys_kos, complements_kos])
   except:
      url_ko_map_colored = "N/A"

   return url_ko_map_colored


def export_pathway_complementarities(config, bins_kos_df):
    """
    Function to get all the KEGG pathway complements among a set of bins

    Input:
    config: output of the config class
    bins_kos_df: dictionary with bin id as a key and the KOs found in the bin as the value

    Returns:
    {beneficiary_bin: {donor_bin_A: {module_a: [], module_b: [],.. }}}
    """
    ko_terms_per_module_definition = config.ko_terms_per_module_definition
    modules_definitions_json_map = config.modules_definitions_json_map
    kegg_maps = config.kegg_modules_to_maps
    maps = open(kegg_maps, "r")
    module_to_map = {}
    for line in maps:
        module, mmap = line.split("\t")
        module_to_map[module[:-1]] = mmap[1:-1]

    # Step 1: Keep track of the KOs related to a module present on each bin
    d = pd.read_csv(ko_terms_per_module_definition, sep="\t")
    d.columns =["module_id","ko_term"]
    d.loc[:, 'presence'] = 1
    definitions_df = d.pivot_table(index='ko_term', columns='module_id', values='presence', fill_value=0)
    ind = definitions_df.index.str.replace('ko:', '')
    definitions_df.index = ind

    bin_kos_per_module = {}
    # Iterate over each column in the second dataframe
    for bin_id in bins_kos_df.columns:
        bin_kos_per_module[bin_id] = {}  # Initialize inner dictionary for each bin
        for module, definition_ko_terms in definitions_df.items():
            # Get KOs of the module definition
            definition_ko_terms = definition_ko_terms[definition_ko_terms != 0]
            # Get KOs present on the bin
            bins_kos = bins_kos_df[bins_kos_df[bin_id]==1][bin_id]
            # Get intersection and add the module: kos_present to the dict
            bin_kos_per_module[bin_id][module] = bins_kos.index.intersection(definition_ko_terms.index).tolist()

    # Step 2: list alternatives for a bin's modules to be completed
    mo_map = json.load(open(modules_definitions_json_map))
    structurals = ["md:M00144","md:M00149","md:M00151",
                   "md:M00152","md:M00154","md:M00155",
                   "md:M00153", "md:M00156", "md:M00158",
                   "md:M00160"]
    # Iterate through bins
    bins_alternatives = {}
    for bin_id in bin_kos_per_module:
        complete_modules = set()
        alternatives_to_gap = {}
        for module, kos_on_its_own in bin_kos_per_module[bin_id].items():
            if module in structurals:
                continue
            # Get KOs related to the module under study that are present on the beneficiary's genome
            list_of_kos_present = set(kos_on_its_own)
            definition_under_study = mo_map[module]['steps']
            definition_under_study_proc = [term if isinstance(term, list) else [term] for term in definition_under_study.values()]
            potential_compl_paths = [list(tup) for tup in itertools.product(*definition_under_study_proc)]
            flat_potent_compl_paths = [flatten(path) for path in potential_compl_paths]
            for path in flat_potent_compl_paths:
                check = all(item in list_of_kos_present for item in path)
                if check:
                    if module not in complete_modules:
                        complete_modules.add(module)
                else:
                    gaps = set(x for x in set(path) if x not in set(list_of_kos_present))
                    if module not in alternatives_to_gap:
                        alternatives_to_gap[module] = {}
                        alternatives_to_gap[module][str(path)] = gaps
                    else:
                        alternatives_to_gap[module][str(path)] = gaps
        # Remove complete modules for the alternatived dict
        for key in complete_modules:
            if key in alternatives_to_gap:
                del alternatives_to_gap[key]
        # Get shortert alternative for each
        for module, path_gaps in alternatives_to_gap.items():
            tmp = tmp2 = alternatives_to_gap[module].copy()
            min_val = min([len(path_gaps[ele]) for ele in path_gaps])
            values = list(tmp2.values())
            shortest_alternatives = [list(tmp2.keys())[values.index(s)]
                                     for s in values
                                     if not any(s.issuperset(i) and len(s) > len(i) for i in values)
                                    ]
            for path, gaps in alternatives_to_gap[module].items():
                if len(gaps) > min_val + 1 or path not in shortest_alternatives:
                    del tmp[path]
            alternatives_to_gap[module] = tmp
        bins_alternatives[bin_id] = alternatives_to_gap

    # Step 3: extract potential complementarities from other bins
    complements = {}
    for beneficiary_bin_id, all_bin_module_alternatives in bins_alternatives.items():
        complements[beneficiary_bin_id] = {}
        for donor_bin_id in bin_kos_per_module:
            complements[beneficiary_bin_id][donor_bin_id] = []
            for module, alts in all_bin_module_alternatives.items():
                donors_kos_relativ_to_module = bin_kos_per_module[donor_bin_id][module]
                for alternative, missing_kos_for_alternative in alts.items():
                    is_subset = set(missing_kos_for_alternative).issubset(set(donors_kos_relativ_to_module))
                    if is_subset:
                        alternative = ast.literal_eval(alternative)
                        beneficiarys_kos_for_alt = set(alternative) - set(missing_kos_for_alternative)
                        try:
                            module_map = module_to_map[module]
                            url = build_kegg_url(module_map, alternative, list(beneficiarys_kos_for_alt))
                        except:
                            url = ""
                        pot_compl = [module,
                                     missing_kos_for_alternative,
                                     alternative,
                                     url
                                     ]
                        complements[beneficiary_bin_id][donor_bin_id].append(pot_compl)
    return bin_kos_per_module, bins_alternatives, complements


def split_list(input_list, chunk_size):
    """
    Split a list to sublists of a size.
    """
    return [input_list[i:i + chunk_size] for i in range(0, len(input_list), chunk_size)]


def run_until_done(command):
    """
    Function to run recursively a command until
    """
    if os.system(command) == 0:
        return 1
    else:
        time.sleep(random.randint(2, 10))
        run_until_done(command)


class build_genres():
    """
    Class to first RAST annotate and the build Genome Scale Metabolic Reconstructions
    using modelseedpy
    """
    def __init__(self, config):
        self.config = config

    def rast_annotate_genomes(self):
        """
        Pool for running rast annotations for a list of bins
        """
        if  len(self.config.bin_filenames) > self.config.threads:
            chunk_size = math.floor( len(self.config.bin_filenames) / self.config.threads  )
            bin_chunks = split_list(self.config.bin_filenames, chunk_size)
        else:
            chunk_size = len(self.config.bin_filenames)
            bin_chunks = [self.config.bin_filenames]

        self.chunk_size = chunk_size

        counter = 0
        for chunk in bin_chunks:
            number_processes = 10
            pool = multiprocessing.Pool(number_processes)
            pool.map(self.rast_annotate_a_genome, chunk)
            pool.close()
            pool.join()
            counter += chunk_size
            print(
                "We now have annotated",
                str(counter),
                "genomes out of the",
                str(len(self.config.bin_filenames))
            )


    def rast_annotate_a_genome(self, bin_filename):
        """
        RAST annotate a user's bin
        """
        bin_file = os.path.join(self.config.bins_path, bin_filename)
        name, _ = os.path.splitext(bin_filename)
        genome = name
        gto_filename = os.path.join(self.config.reconstructions, "".join([name, ".gto"]))

        # rast_create_genome_command
        rast_create_genome_command = " ".join([
            "rast-create-genome", "--scientific-name", name,
            "--genetic-code", "11", "--domain", "Bacteria",
            "--contigs", bin_file,
            "--genome-id", genome, ">", gto_filename
            ])
        if not os.path.exists(gto_filename):
            print("rast_create_genome_command:", rast_create_genome_command)
            run_until_done(rast_create_genome_command)

        # rast-process-genome
        gto_filename_2 = "".join([gto_filename, "_2"])
        rast_process_genome_command = " ".join([
            "rast-process-genome", "<", gto_filename, ">", gto_filename_2
            ])
        if not os.path.exists(gto_filename_2):
            print("\nrast_process_genome_command: ", rast_process_genome_command)
            run_until_done(rast_process_genome_command)

        # rast-export-genome protein_fasta
        faa_filename = "".join([gto_filename[:-4], ".faa"])
        rast_export_genome_command = " ".join([
            "rast-export-genome",
            "protein_fasta",
            "<", gto_filename_2, ">",
            faa_filename
            ])
        if not os.path.exists(faa_filename):
            print("\nrast_export_genome_command: ", rast_export_genome_command)
            run_until_done(rast_export_genome_command)

        # rast-export-genome seed_dir
        seed_dir = "".join([gto_filename[:-4], ".seed_dir.tar.gz"])
        rast_export_genome_dir_command = " ".join([
            "rast-export-genome", "seed_dir",
            "<", gto_filename_2, ">",
            seed_dir
            ])
        if not os.path.exists(seed_dir):
            print("\nrast_export_genome_dir_command:", rast_export_genome_dir_command)
            run_until_done(rast_export_genome_dir_command)

        return


    def modelseed_reconstructions(self):
        """
        Pool for running GENREs reconstruction using the final .faa files from the
        rast_annotate_genomes() function
        [NOTE] Not to be used for now as the RAST server seems not that stable to have several queries..
        """
        faa_files = [file
                     for file in os.listdir(self.config.reconstructions)
                     if file.endswith(".faa")
                    ]
        # if len(faa_files) > self.config.threads:
        #     faa_chunks = bin_chunks = split_list(faa_files, self.chunk_size)
        # else:
        #     faa_chunks = [faa_files]

        # counter = 0
        # for chunk in faa_chunks:
        #     pool = multiprocessing.Pool(self.config.threads)
        #     pool.map(self.reconstruct_a_model, chunk)
        #     pool.close()
        #     pool.join()
        #     counter += self.chunk_size
        #     print(
        #         "We now have reconstructed",
        #         str(counter),
        #         "GENREs out of the",
        #         str(len(faa_files))
        #     )
        for faa_file in faa_files:
            print("faa file to be reconstructed:", faa_file)
            self.reconstruct_a_model(faa_file)

    def reconstruct_a_model(self, annotation_faa):
        """
        Build a Genome Scale reconstruction using ModelSEEDpy and the PATRIC annotations
        """
        # ModelSEED using the default b.f and complete medium, gapfill with the default algo
        model_id, _ = os.path.splitext(annotation_faa)
        annotation_faa_path = os.path.join(self.config.reconstructions, annotation_faa)
        model_filename =  os.path.join(self.config.genres, "".join([model_id, ".xml"]))

        if os.path.exists(model_filename):
            return 1

        msgenome = MSGenome.from_fasta(annotation_faa_path, split=' ')
        print("msgenome ready for build_metabolic_model() ")
        model = self.recursive_build(model_id, msgenome)

        cobra.io.write_sbml_model(cobra_model = model, filename = model_filename)

        return 1

    def recursive_build(self, model_id, msgenome):
        """
        RAST server usually has issues that lead to fail attempts.
        This recursive function allows the build_metabolic_model() to be performed until the server responses.
        """
        try:
            model = MSBuilder.build_metabolic_model(model_id = model_id,
                                                    genome   = msgenome,
                                                    index    = "0",
                                                    # classic_biomass = True,
                                                    gapfill_model   = self.config.gapfill_model,
                                                    gapfill_media   = None,
                                                    annotate_with_rast = True,
                                                    allow_all_non_grp_reactions = True
                                                    )
            return model
        except:
            time.sleep(random.randint(5, 15))
            return self.recursive_build(model_id, msgenome)


def get_seed_sets(config):
    """
    """
    all_genres_files = [os.path.join(config.reconstructions, file) for file in os.listdir(config.reconstructions)]
    genre_files = [file for file in all_genres_files if file.startswith(".xml")]
    for file in genre_files:
        dest_path = os.path.join(config.genres, os.path.basename(file))
        shutil.move(file, dest_path)

    phylomint_params = ["./PhyloMInt",
                        "-d", config.genres,
                        "--outdir", config.seeds,
                        "-o", "phylomint_scores.tsv",
                        "--dics", "True",
                        "--threads", str(config.threads)
                        ]
    phylomint_cmd = " ".join(phylomint_params)
    print(phylomint_cmd)
    os.system(phylomint_cmd)

