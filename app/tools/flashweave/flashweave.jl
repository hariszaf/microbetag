# microbetag is working as a set of workflows that consists of tools
# This script runs FlashWeave to get a co-occurrence network based on 
# an OTU or ASV table (mandatory) and a metadata table (optional)
# author: Haris Zafeiropoulos
# date  : 2022-01

# Essential parameters:
#     n_obs_min - don't compute associations between variables having less reliable samples (i.e. non-zero if heterogeneous=true) than this number.
#                   -1: automatically choose a threshold.
#     track_rejections - store for each discarded edge, which variable set lead to its exclusion (can be memory intense for large networks)
#     transposed - if true, rows of data are variables and columns are samples



using FlashWeave # this has some pre-compilation delay the first time it's called, subsequent imports are fast


if length(ARGS)>1
   data_path = ARGS[1]
   metadata = ARGS[2]
   netw_results = learn_network(
                                 data_path, 
                                 meta_data_path, 
                                 sensitive = true, 
                                 heterogeneous = false
                              )


else
   data_path = ARGS[1]
   netw_results = learn_network(
                                 data_path,
                                 n_obs_min = 5,  
                                 sensitive = true, 
                                 heterogeneous = false,
                                 max_k = 1,
                                 transposed = true
                              )

end



G = graph(netw_results) # weighted graph object representing interactions + weights, 
                        # to be used with the JuliaGraphs ecosystem (https://github.com/JuliaGraphs)

save_network("/mnt/network_output.edgelist", netw_results)
save_network("/mnt/network_detailed_output.edgelist", netw_results, detailed=true)

