#!/usr/bin/env python
# Copyright Rogan Carr, Roie Levy, and Elhanan Borenstein 2013
#  Contact: rogan@uw.edu
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# SupportingNodes == nodesAll -> BSS; == nonSeeds -> metabolic_complementarity
def compute_single_interaction_score(SeedGroups,SupportingNodes):
    
    score=0
    SupportedNodes=[]
    
    for group in SeedGroups:
        found=0
        for seed in SeedGroups[group]:
            if seed in SupportingNodes:
                if not found:
                    score+=1
                    found=1
                SupportedNodes.append(seed)
                #break
            
    score /= float(len(SeedGroups))

    return score, SupportedNodes

def compute_all_interaction_scores(SeedGroupOne,AllNodesOne,NonSeedsOne,SeedGroupTwo,AllNodesTwo,NonSeedsTwo):
    
    # Biosynthetic support score is seeds vs. All nodes
    BSSOneTwo=compute_single_interaction_score(SeedGroupOne,AllNodesTwo)
    BSSTwoOne=compute_single_interaction_score(SeedGroupTwo,AllNodesOne)
    
    # Metabolic complementarity score is seeds vs. non-seed nodes
    MCOneTwo=compute_single_interaction_score(SeedGroupOne,NonSeedsTwo)
    MCTwoOne=compute_single_interaction_score(SeedGroupTwo,NonSeedsOne)
    
    return BSSOneTwo, BSSTwoOne, MCOneTwo, MCTwoOne
