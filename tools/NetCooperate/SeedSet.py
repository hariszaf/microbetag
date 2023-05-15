#!/usr/bin/env python
# Copyright Rogan Carr and Elhanan Borenstein 2013
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

# Calculate the seed set

# Do we really need the node order?
# I think it's actually used later to define which nodes are in the search
# Yeah, it's arbitrary
def calculate_seeds(graph,onlyGiant=None,minComponentSize=None,seedThreshold=0):


    # Create our datastructures

    #  Construct a reverse graph
    reverse_graph={}

    #  Calculate the node order
    nodes=set()
    nodeOrder=[];
    
    for source_node in graph:

        # Add it to the list of nodes
        if not source_node in nodes:
            nodes.add(source_node)
            nodeOrder.append(source_node)

        for target_node in graph[source_node]:
            if not reverse_graph.has_key(target_node):
                reverse_graph[target_node]=[source_node]
            else:
                reverse_graph[target_node].append(source_node)

            # Add it to the list of nodes
            if not target_node in nodes:
                nodes.add(target_node)
                nodeOrder.append(target_node)
    
    # We now return all the nodes that were not pruned
    #del nodes
    

    Pruned=set()
    # Trim (all or some) Small components
    if onlyGiant or minComponentSize:
        
        # Use DFS to determine the connected components
        nComponents,nodeToComponent = depth_first_search(graph,nodeOrder,reverse_graph)[1:3]
        
        # Calculate the size of each component
        componentSize = [0]*nComponents        

        for node in nodeOrder:
            componentSize[nodeToComponent[node]]+=1

        # If we only keep the giant component
        #  Then set the minimum size to the giant c's size
        if onlyGiant:
            minComponentSize=0
            for c in xrange(nComponents):
                if componentSize[c] > minComponentSize:
                    minComponentSize=componentSize[c]
        
        # Delete all nodes not in the minimal component

        #  Identify the components to trim
        componentsToTrim = set()
        for c in xrange(nComponents):
            if componentSize[c] < minComponentSize:
                componentsToTrim.add(c)
        
        #  Trim the nodes
        #    Note the .keys() function is necessary because we are deleting keys!
        #    Also, we add to Pruned, remove frome nodeOrder
        for node in graph.keys():
            if nodeToComponent[node] in componentsToTrim:
                del graph[node]

                if node in reverse_graph:
                    del reverse_graph[node]
                Pruned.add(node)

                if node in nodeOrder:
                    nodeOrder.remove(node)
        
        #  Capture the only-pointed-to nodes
        #    Note the .keys() function is necessary because we are deleting keys!
        for node in reverse_graph.keys():
            if nodeToComponent[node] in componentsToTrim:
                del reverse_graph[node]
                Pruned.add(node)
                if node in nodeOrder:
                    nodeOrder.remove(node)
        
        # We return all the nodes (it's the convention) 
        #  #Trim the master list of nodes
        #for node in Pruned:
        #    nodes.discard(node);

    # Identify the Strongly-Connected Components (SCC)
    finishOrder,nComponents = depth_first_search(graph,nodeOrder)[0:2]
    nSCC,nodeToSCC = depth_first_search(reverse_graph,finishOrder)[1:3]
    
    # Identify seed sets

    #  Calculate the size of each scc
    sccSize = [0]*nSCC
    for node in nodeOrder:
        sccSize[nodeToSCC[node]]+=1
    
    # Determine which SCCs are sources
    #   Initialize to all SCCs
    SeedGroups={}
    for i in xrange(nSCC):
        SeedGroups[i]=[]
    #   Remove SCCs that are targets of others
    for source_node in graph:
        for target_node in graph[source_node]:
            if nodeToSCC[target_node] != nodeToSCC[source_node]:
                if nodeToSCC[target_node] in SeedGroups:
                    del SeedGroups[nodeToSCC[target_node]]

    # Assemble the seeds
    Seeds={}
    nonSeeds={}
    for node in nodeOrder:
        seed_rank=0
        if nodeToSCC[node] in SeedGroups:
            seed_rank = 1 / float(sccSize[nodeToSCC[node]])
            
        if seed_rank > seedThreshold:
            Seeds[node]=seed_rank
            SeedGroups[nodeToSCC[node]].append(node)
        else:
            # This can be nonzero
            nonSeeds[node]=seed_rank

    # Remove seed groups that are defunct due to their size
    #    Note the .keys() function is necessary because we are deleting keys!
    for group in SeedGroups.keys():
        if len(SeedGroups[group]) < sccSize[group]:
            if len(SeedGroups[group]) == 0 and 1/float(sccSize[group]) <= seedThreshold:
                print "Removing Seed Group", group, "of size", sccSize[group]
                del SeedGroups[group]
            else:
                print "Serious problem in SeedGroup", group, ": Expected", sccSize[group], "and found", len(SeedGroups[group]),".";
            
    return [Seeds, SeedGroups, nonSeeds, Pruned, nodes]
        
# Particular to Our Strategy
def depth_first_search(graph,nodeOrder,reverse_graph=[]):
    
    # To return
    component=-1
    nodeToComponent={}
    finishOrder=[]
    
    # Our search space
    stack=[]

    # Working
    visited={}
    
    # Search every node
    for node in reversed(nodeOrder):
        if node in visited:
            continue
        
        # We must be on a new component if we haven't been here
        component+=1

        # Put it on our search space
        stack.append(node)

        # ## to here (from below)
        #RC swapped this out
        
        # # Assign the compoment
        # nodeToComponent[node]=component
                
        # # Mark it as visited
        # visited[node]=1
        
        # # Mark the neighbors
        # dfs_get_neighbors(graph,visited,node,stack)
        # # Get the neighobrs for the revers graph if it's defined
        # if reverse_graph:
        #     dfs_get_neighbors(reverse_graph,visited,node,stack)
        
        # # I think I can roll this up into the loop from here ^^

        # We visit each node twice
        while stack:
            vertex = stack.pop()

            # If it's visited, roll it up
            if vertex in visited:
                if visited[vertex]<2:
                    finishOrder.append(vertex)
                    visited[vertex]=2
            else:
                # Assign the component
                nodeToComponent[vertex] = component
                # Put it on our search space
                stack.append(vertex)
                
                # Mark it as visited
                visited[vertex]=1

                dfs_get_neighbors(graph,visited,vertex,stack)
                # Get the neighobrs for the revers graph if it's defined
                if reverse_graph:
                    dfs_get_neighbors(reverse_graph,visited,vertex,stack)

    # convert to a 1-based list
    component+=1
    return [finishOrder, component, nodeToComponent]

def dfs_get_neighbors(graph,visited,node,stack):
    if node in graph:
        for next_node in graph[node]:
            if next_node not in visited:
                stack.append(next_node)

