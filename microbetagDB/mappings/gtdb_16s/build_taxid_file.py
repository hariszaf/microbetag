import sys 

ranks = open("gtdb_full_taxonomies.txt", "r")
ranks = ranks.readlines()
ranks = [line.strip().split(";") for line in ranks]

taxa = [["domain", "phylum", "class", "order", "family", "genus", "species"],
        ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]]

count = 1
groups = ["Root"]
index = -1
level = 0
rank = "rootrank"


myDict = {}
myDict[0] = {"group": "Root", "parent": -1, "rank_index":0, "level": "rootrank"}
myDict[1] = {"group": "Bacteria", "parent": 0, "rank_index":1, "level": "domain"}
myDict[2] = {"group": "Archaea", "parent": 0, "rank_index":1, "level": "domain"}

combos = {("Root","rootrank"):0, ("Bacteria", "domain"): 1, ("Archaea", "domain"): 2}

import sys

for i in range(len(ranks)):

    if i%1000 == 0:
        print(i, "out of", str(len(ranks)))

    for j in range(len(ranks[i])):

        index = len(myDict)
        group = ranks[i][j][3:]

        level_short = ranks[i][j][:3]
        rank_index = taxa[1].index(level_short)
        level = taxa[0][rank_index]

        if j >= 1:

            parent_group = ranks[i][j-1][3:]
            parent_level_short = taxa[1][rank_index -1]
            parent_level = taxa[0][rank_index -1]
            parent_rank_index = taxa[1].index(parent_level_short)

            if ((parent_group, parent_level)) in combos.keys():

                parent = combos[(parent_group, parent_level)]
        
                cand = {"group": group, "parent": parent, "rank_index": rank_index+1, "level": level}

                if cand not in myDict.values():
                    myDict[index] = cand
                    combos.update({(group, level): index})

with open("gtdb_taxid.txt", "w") as f:
    for k,v in myDict.items():
        t = "".join([str(k), "*", "*".join( [str(x) for x in list(v.values())])])
        f.write(t + "\n")




