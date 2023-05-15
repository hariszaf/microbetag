
#gtdb_core_dic = {}
#ncbi_core_dic = {}
#ncbi_unf_dic  = {}
#silva_ssu_dic = {}

gtdb_core_file = open("gtdb_core", "r")
#species = open("species2ncbiId.tsv","r")
output = open("alternative", "w")

#for line in gtdb_core_file: 
#   line = line.split("\t")
#   gtdb_core_dic[line[4]] = line[1]
#   ncbi_core_dic[line[3]] = line[1]
#   ncbi_unf_dic[line[6][:-1]] = line[1]
#   silva_ssu_dic[line[5]] = line[1]
#
#
#for x in species:
#   name = x.split("\t")[0]
#   ncbiId = x.split("\t")[1][:-1]
#   if name in gtdb_core_dic.keys():
#      output.write(name + "\t" + ncbiId + "\t" + gtdb_core_dic[name] + "\n")
#   elif name in ncbi_core_dic.keys():
#      output.write(name + "\t" + ncbiId + "\t" + ncbi_core_dic[name] + "\n")
#   elif name in ncbi_unf_dic.keys():
#      output.write(name + "\t" + ncbiId + "\t" + ncbi_unf_dic[name] + "\n")
#   elif name in silva_ssu_dic.keys():
#      output.write(name + "\t" + ncbiId + "\t" + silva_ssu_dic[name] + "\n")
#   else:
#      print("No GTDB ref genome for species: ", name)
#


two_cols = open("species2ncbiId.tsv", "r")
two_cols = two_cols.readlines()
three_cols = open("species2ncbiId2gtdbGenome.tsv","r")
three_cols = three_cols.readlines()
two_cols_dic = {}

for entry in two_cols:
   id = entry.split("\t")[1][:-1]
   name = entry.split("\t")[0]
   if id not in two_cols_dic:
      two_cols_dic[id] = [name]
   else:
       two_cols_dic[id].append(name)

three_cols_ids = {}
three_cols_names = []
for entry in three_cols:
   id = entry.split("\t")[1]
   accession = entry.split("\t")[2][:-1]
   name = entry.split("\t")[0]
   if id not in three_cols_ids:
      three_cols_ids[id] = accession
   three_cols_names.append(name)

output = open("more.tsv", "w")
for ncbiId, names in two_cols_dic.items():
   if ncbiId in three_cols_ids:
      for name in names: 
         if name not in three_cols_names:
            output.write(name + "\t" + ncbiId + "\t" + three_cols_ids[ncbiId] + "\n")


gtdb_core = open("gtdb_core","r")
qiime = open("dada2species.tsv", "r")
output = open("dada2ncbi2accession.tsv", "w")
ncbiId_accession = {}
silva_accession = {}
silva_unf_accession = {}

for line in gtdb_core:
   line = line.split("\t")
   accession = line[1]
   ncbiId = line[2]
   silva = line[3]
   silva_unf = line[4][:-1]
   ncbiId_accession[ncbiId] = accession
   silva_accession[silva] = accession
   silva_unf_accession[silva_unf] = accession

for line in qiime:
   line = line.split("\t")
   ncbiId = line[1][:-1]
   name = line[0]
   if ncbiId in ncbiId_accession:
      output.write( name + "\t" + ncbiId + "\t" + ncbiId_accession[ncbiId] + "\n" )
   elif name in silva_accession:
      output.write( name + "\t" + ncbiId + "\t" + silva_accession[name] + "\n" )
   elif name in silva_unf_accession: 
      output.write( name + "\t" + ncbiId + "\t" + silva_unf_accession[name] + "\n" )      
   else:
      output.write( name + "\t" + ncbiId + "\tnan\n")


overall_file = open("species2ncbiId2accession.tsv", "r")
triplets = {}
taxid_accession = {}
for line in overall_file:
   line = line.split("\t")
   name = line[0]
   taxid = line[1]
   gca = line[2][:-1]
   if name not in triplets:
      triplets[name] = [taxid]
   else:
      triplets[name].append(taxid)
   if taxid not in taxid_accession:
      taxid_accession[taxid] = [gca]
   else:
      taxid_accession[taxid].append(gca)
print(triplets["Ruminococcus flavefaciens"])


for name, ids in triplets.items():
   ids_counts = {}
   for id in ids:
      if id not in ids_counts:
         ids_counts[id] = 1
      else: 
         ids_counts[id] += 1
   if len(ids_counts) > 1:
      x = list(ids_counts.values())
      if([x[0]]*len(x) == x):
         # keep them all 
         print("domage: ", name)
         print(ids_counts)
         print("~~~~")
      else:
         key_id = max(ids_counts, key=ids_counts.get)
         print("bien sure: ", name)
         print(key_id)
         print("~~~~")


outut = open("mplasdasdasd.tsv", "w")





