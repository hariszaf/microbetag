# author: Haris Zafeiropoulos

import json
import sys
import os
import subprocess
from argparse import ArgumentParser
from pathlib import Path

# Libraries
from libsbml import readSBML

# Local modules
from lib import BuildGraphNetX, CalculateIndexes

subf = str(sys.argv[1])

dirr_path="/staging/leuven/stg_00106/haris/microbetag/RECONSTRUCTIONS/all_recons/subfolder_" + subf

# get all XML files in directory
sbml_files = [ f for f in os.listdir(dirr_path) if f.endswith('.xml') ]
temp_set = set()
SeedSetDic = dict()
nonSeedSetDic = dict()
ConfidenceDic = dict()
count = 0
total = len(sbml_files)
for sbml in sbml_files:
    count += 1
    sbml_path = ('%s/%s' % (dirr_path, sbml))
    sbml_base = sbml.rstrip('.xml')
    print(f'Calculating Seeds: {sbml_base} - {count}/{total}')
    if sbml_base not in temp_set:
        # calculate SeedSets
        DG_sbml = BuildGraphNetX.buildDG(sbml_path)
        SeedSetConfidence, SeedSet, nonSeedSet = BuildGraphNetX.getSeedSet(DG_sbml, maxComponentSize=5) # args.maxcc default value
        # Save dicts
        SeedSetDic[sbml_base] = SeedSet
        nonSeedSetDic[sbml_base] = nonSeedSet
        ConfidenceDic[sbml_base] = SeedSetConfidence
        temp_set.add(sbml_base)

with open("subfolder_" + subf + "_ConfDic.json", "w") as json_file:
    json.dump(ConfidenceDic, json_file)
with open("subfolder_" + subf + "_nonSeedSetDic.json", "w") as json_file:
    json.dump(nonSeedSetDic, json_file)
