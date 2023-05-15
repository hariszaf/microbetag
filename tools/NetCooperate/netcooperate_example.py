#!/usr/bin/env python

from __future__ import print_function
import sys
import SeedSet
import NetCooperate

if len(sys.argv) != 3:
    print(sys.argv[0], 'graph1 graph2');
    exit();

# Our graph reading script
#  Reads a tab-separated file of edges
def readGraph(graphFile):
    
    graph={}
    for line in open(graphFile):
        
        vals = line.strip().split("\t");
        
        if len(vals[0])==0 or vals[0][0] == '#':
            continue     
        elif len(vals) != 2:
            print("Bad line: " + line)
            continue
        
        graph.setdefault(vals[0],[]).append(vals[1])
        
    return graph

# Read in the files
graphFileOne=sys.argv[1];
graphFileTwo=sys.argv[2];

# Get our graphs
graphOne=readGraph(graphFileOne);
graphTwo=readGraph(graphFileTwo);

# Choose what we find important
onlyGiant=False
minComponentSize=0

# Calculate the seed sets
SeedsOne, SeedGroupsOne, nonSeedsOne, PrunedOne, NodesOne = SeedSet.calculate_seeds(graphOne,onlyGiant,minComponentSize)
SeedsTwo, SeedGroupsTwo, nonSeedsTwo, PrunedTwo, NodesTwo = SeedSet.calculate_seeds(graphTwo,onlyGiant,minComponentSize)

# Calculate the biosynthetic support
#  The biosynthetic support compares the SeedGroups to All nodes in the other network

bssOneOne, bssOneOneSupport = NetCooperate.compute_single_interaction_score(SeedGroupsOne,NodesOne)
bssOneTwo, bssOneTwoSupport = NetCooperate.compute_single_interaction_score(SeedGroupsOne,NodesTwo)
bssTwoOne, bssTwoOneSupport = NetCooperate.compute_single_interaction_score(SeedGroupsTwo,NodesOne)
bssTwoTwo, bssTwoTwoSupport = NetCooperate.compute_single_interaction_score(SeedGroupsTwo,NodesTwo)

# Calculate the metabolic complementarity
#  The metabolic complementarity compares the SeedGroups to Non-Seed nodes in the other network

mcOneOne, mcOneOneSupport = NetCooperate.compute_single_interaction_score(SeedGroupsOne,nonSeedsOne)
mcOneTwo, mcOneTwoSupport = NetCooperate.compute_single_interaction_score(SeedGroupsOne,nonSeedsTwo)
mcTwoOne, mcTwoOneSupport = NetCooperate.compute_single_interaction_score(SeedGroupsTwo,nonSeedsOne)
mcTwoTwo, mcTwoTwoSupport = NetCooperate.compute_single_interaction_score(SeedGroupsTwo,nonSeedsTwo)

print()
print("Example one: computing single interaction scores individually & reporting supported compounds.\n")
print("biosynthetic support:\n")
print("%s being supported by %s:\t%g\t%s" % (graphFileOne, graphFileOne, bssOneOne, ','.join(bssOneOneSupport)))
print("%s being supported by %s:\t%g\t%s" % (graphFileOne, graphFileTwo, bssOneTwo, ','.join(bssOneTwoSupport)))
print("%s being supported by %s:\t%g\t%s" % (graphFileTwo, graphFileOne, bssTwoOne, ','.join(bssTwoOneSupport)))
print("%s being supported by %s:\t%g\t%s" % (graphFileTwo, graphFileTwo, bssTwoTwo, ','.join(bssTwoTwoSupport)))
print()
print("metabolic complementarity:\n")
print("%s being complemented by %s:\t%g\t%s" % (graphFileOne, graphFileOne, mcOneOne, ','.join(mcOneOneSupport)))
print("%s being complemented by %s:\t%g\t%s" % (graphFileOne, graphFileTwo, mcOneTwo, ','.join(mcOneTwoSupport)))
print("%s being complemented by %s:\t%g\t%s" % (graphFileTwo, graphFileOne, mcTwoOne, ','.join(mcTwoOneSupport)))
print("%s being complemented by %s:\t%g\t%s" % (graphFileTwo, graphFileTwo, mcTwoTwo, ','.join(mcTwoTwoSupport)))

### Or use this all-in-one function call

# This returns a two-tuple for each argument. tupple[0] is the support score, tupple[1] is the list of supported seeds
BSSOneOnTwo, BSSTwoOnOne, MCOneOnTwo, MCTwoOnOne = NetCooperate.compute_all_interaction_scores(SeedGroupsOne,NodesOne,nonSeedsOne,SeedGroupsTwo,NodesTwo,nonSeedsTwo)

print("\n---------------------------------------------------\n")
print("Example two: computing all interaction scores in one call but not reporting supported compounds.\n")
print("biosynthetic support:\n")
print("%s being supported by %s:\t%g" % (graphFileOne, graphFileOne, BSSOneOnTwo[0]))
print("%s being supported by %s:\t%g" % (graphFileOne, graphFileTwo, BSSTwoOnOne[0]))
print("")
print("metabolic complementarity:\n")
print("%s being complemented by %s:\t%g" % (graphFileTwo, graphFileOne,  MCOneOnTwo[0]))
print("%s being complemented by %s:\t%g" % (graphFileTwo, graphFileTwo,  MCTwoOnOne[0]))

