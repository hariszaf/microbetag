#!/usr/bin/env python
import yaml 
import os
import pandas as pd


class Config:
    def __init__(self, io_path, conf):
        self.io_path = io_path
        self.abundance_table = os.path.join(io_path, conf["abundance_table_file"]["path"])
        self.build_network = conf["build_network"]["value"]
        # User's group and user id
        self.user_id = int(os.getenv('USER_ID'))
        self.group_id = int(os.getenv('HOST_GID'))
        # Output dir
        self.output_dir = os.path.join(self.io_path, conf["output_directory"]["path"])
        # Metadata file for FlashWeave
        if conf["metadata_file"]["path"]:
            self.metadata = "true"
            self.metadata_file = os.path.join(io_path, conf["metadata_file"]["path"])
        else:
            self.metadata = "false"
        # Flashweave arguments
        if conf["flashweave_sensitive"]["value"]:
            self.sensitive = "true"
        else:
            self.sensitive = "false"
        if conf["flashweave_heterogeneous"]["value"]:
            self.heterogeneous = "true"
        else:
            self.heterogeneous = "false"

# Vars
io_path = "/media"
main_path = "/pre_microbetag"
flashweave_input_file = "abundance_table_flashweave.tsv"
sequence_file_for_taxid = "seqs_tmp.fa"
tax_assignmets_file = "assignments.tsv"
julia_path = "/opt/julia-1.7.1/bin/julia/julia"
rscript_path = "/usr/local/bin/Rscript"
classify_Rscript = os.path.join(main_path, "classify.R")
flashweave_script = os.path.join(main_path, "flashweave.jl")

config_file = os.path.join(io_path, "config.yml")
with open(config_file, 'r') as yaml_file:
    config = Config(io_path, yaml.safe_load(yaml_file))

try:
    os.mkdir(config.output_dir)
except FileExistsError as e:
    warning_message = f"Warning: The output directory '{config.output_dir}' already exists. Its content will be overwritten with the new - {e}"
    print(warning_message)

if config.abundance_table:

    """check if abundance table there and a sequence column included"""
    if config.abundance_table is None:
        raise ValueError("You need to provide an abundance table including a column with the ASV/OTU sequence.")

    """get taxonomy assignments using DECIPHER and 16S GTDB ref seqs"""
    abundance_table_data = pd.read_csv(config.abundance_table, sep=None,  engine='python') 
    column_names = list(abundance_table_data.columns)
    seqid = column_names[0]
    seq = column_names[-1]
    seqs_df = abundance_table_data[ [seqid, seq] ]
    with open(sequence_file_for_taxid, "w") as f:
        for pair in seqs_df.to_dict(orient="records"):
            line_to_write = "".join([">", pair[seqid], "\n", pair[seq], "\n"])
            f.write(line_to_write)

    # fix Rscript path 
    classify_command = " ".join([rscript_path, classify_Rscript])
    if os.system(classify_command) != 0:
        raise SystemError("Taxonomy assignment using IDTAXA of the DECIPHER package and the 16S GTDB sequences as reference failed.\
                          Please check R dependencies are ok and that you have the trained gtdb_16S.RData file.")

    assignments = pd.read_csv(tax_assignmets_file, sep="\t")
    merged_df = pd.merge(abundance_table_data, assignments, left_on=seqid, right_on='seqid', how='left')

    # Drop the redundant 'taxid' column
    merged_df = merged_df.drop('seqid', axis=1)

    gtdb_df = merged_df.drop(seq, axis=1)
    tax_annotated_abundance_table_file = os.path.join(config.output_dir, "GTDB_tax_assigned_abundance_table.tsv")
    gtdb_df.to_csv(tax_annotated_abundance_table_file, sep="\t", index=False)


if config.build_network:

    flashweave_table = abundance_table_data.drop(seq, axis = 1)
    float_col        = flashweave_table.select_dtypes(include=['float64']) 
   
    for col in float_col.columns.values:
        flashweave_table[col] = flashweave_table[col].astype('int64')

    flashweave_table[seqid] = flashweave_table[seqid].astype(str)
    flashweave_table.to_csv(flashweave_input_file, sep='\t', index=False)

    print("A cooccurrence network is about to be built using FlashWeave.")

    if config.metadata == "true":
        flashweave_params = [
            julia_path, flashweave_script, config.output_dir, flashweave_input_file, config.sensitive, config.heterogeneous, config.metadata, config.metadata_file
        ]
    else:
        flashweave_params = [
            julia_path, flashweave_script, config.output_dir, flashweave_input_file, config.sensitive, config.heterogeneous, config.metadata
        ]
    flashweave_command = " ".join(flashweave_params)

    if os.system(flashweave_command) != 0:
        raise SystemError("FlashWeave failed. Check if Julia and FlashWeave is there. If yes, check the abundance table you provide.")

# Change ownership to output folder and its contents
os.chown(config.output_dir, config.user_id, config.group_id)
for filename in os.listdir(config.output_dir):
    file_path = os.path.join(config.output_dir, filename)
    os.chown(file_path, config.user_id, config.group_id)

