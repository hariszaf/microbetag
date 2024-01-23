# author: Haris Zafeiropoulos

import json
import sys
import os
import subprocess
from argparse import ArgumentParser
from pathlib import Path
import multiprocessing
import time

# Libraries
from libsbml import readSBML

# Local modules
from lib import BuildGraphNetX, CalculateIndexes

all_recon_path = "/staging/leuven/stg_00106/haris/microbetag/RECONSTRUCTIONS/all_recons"
parsed_ids = open(all_recon_path + "/completed.txt", "r").readlines()
parsed_ids = [c[:-1] for c in parsed_ids]
all_recon_files = os.listdir(all_recon_path)

ConfidenceDic = {}
nonSeedSetDic = {}
sbml_files = []
for ffile in all_recon_files:
	if ffile.endswith("ConfDic.json"):
		ConfidenceDic.update(json.load(open(ffile, "r")))
	elif ffile.endswith("nonSeedSetDic.json"):
		nonSeedSetDic.update(json.load(open(ffile, "r")))
	elif os.path.isdir(ffile):
		if ffile.startswith("subfolder"):
			sbml_files.extend([ f for f in os.listdir(ffile) if f.endswith('.xml.gz') ])

total = len(sbml_files)

#print("all SBML files:", str(total))
#print("parsed ids: ", str(len(parsed_ids)))

count = 0
pairwise_total = total*total

def apply_scores(A):
	from lib import BuildGraphNetX, CalculateIndexes
	A = A.rstrip('.xml.gz')
	if A in parsed_ids:
		print(A, " already parsed")
		return
	print(A, "not parsed yet!")
	outs = []
	for B in sbml_files:
		B = B.rstrip('.xml.gz')
		MetabolicCompetitionIdxAB = CalculateIndexes.MetabolicCompetitionIdx(ConfidenceDic[A], ConfidenceDic[B])
		MetabolicCooperationIdxAB = CalculateIndexes.MetabolicCooperationIdx(ConfidenceDic[A], ConfidenceDic[B], nonSeedSetDic[B])
		output = (f'{A}\t{B}\t{MetabolicCompetitionIdxAB}\t{MetabolicCooperationIdxAB}\n')
		#out_file.write(output)
		#print(output)
		outs.append(output)
	return outs

"""
with open('phylomint_all_scores.tsv', 'a') as out_file:
	#out_file.write(f'\tA\tB\tCompetition\tComplementarity\n')
	for A in sbml_files:
		count += 1
		A = A.rstrip('.xml.gz')
		print(A)
		if A in parsed_ids:
			print("Already parsed")
			continue
		print("Getting scores for file ", A, " number ", str(count), " out of ", str(total))
		for B in sbml_files:
			# A = A.rstrip('.xml')
			B = B.rstrip('.xml.gz')
			#count += 1
			#print(f'Calculating Indexes: {A} vs {B} - {count}/{pairwise_total}')
			MetabolicCompetitionIdxAB = CalculateIndexes.MetabolicCompetitionIdx(ConfidenceDic[A], ConfidenceDic[B])
			MetabolicCooperationIdxAB = CalculateIndexes.MetabolicCooperationIdx(ConfidenceDic[A], ConfidenceDic[B], nonSeedSetDic[B])
			output = (f'{A}\t{B}\t{MetabolicCompetitionIdxAB}\t{MetabolicCooperationIdxAB}\n')
			out_file.write(output)

"""




sbml_files_hundreds =  [sbml_files[i:i+500] for i in range(0, len(sbml_files), 500)]

counterr = 0

print("Number of iterations to perform: ", str(len(sbml_files_hundreds)))

for i  in range(len(sbml_files_hundreds)): # (len(sbml_files_hundreds)):

	print("\n\n\n  >>>>>>>>>>>>>   hundred: ", str(counterr), "   <<<<<<<<<<<<<<")
	start_time = time.time()
	pool = multiprocessing.Pool(30)
	results = pool.map(apply_scores, sbml_files_hundreds[0])
	print("Write scores for running iteration\n")
	with open('phylomint_all_scores.tsv', 'a') as out_file:
		for species in results:
			# If a species is already parsed, it returns nothing
			if species is not None:
				for line in species:
					# In case anythinge is empty..
					if line is not None:
						out_file.write(line)

	counterr += 1
	pool.close()
	pool.join()

	end_time = time.time()
	execution_time = end_time - start_time
	print(f"Iteration execution time: {execution_time:.6f} seconds")
