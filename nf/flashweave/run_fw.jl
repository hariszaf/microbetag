using ArgParse
using FlashWeave

# Parse arguments
s = ArgParseSettings()
# @add_arg_table s begin
#     "--input"        , arg_type = String,  help    = "Input abundance table"
#     "--metadata"     , arg_type = String,  default = "", help = "Optional metadata file"
#     "--sensitive"    , arg_type = Bool,    default = true
#     "--heterogeneous", arg_type = Bool,    default = false
#     "--max_k"        , arg_type = Int,     default = 3
#     "--n_obs_min"    , arg_type = Int,     default = -1
#     "--alpha"        , arg_type = Float64, default = 0.01
# end
@add_arg_table s begin
    "--input"
        help = "Input abundance table"
        arg_type = String
        required = true

    "--metadata"
        help = "Optional metadata file"
        arg_type = String
        default = ""

    "--sensitive"
        help = "Use sensitive mode"
        arg_type = Bool
        default = true

    "--heterogeneous"
        help = "Heterogeneous network"
        arg_type = Bool
        default = false

    "--max_k"
        help = "Max k neighbors"
        arg_type = Int
        default = 3

    "--n_obs_min"
        help = "Minimum number of observations"
        arg_type = Int
        default = -1

    "--alpha"
        help = "Significance level"
        arg_type = Float64
        default = 0.01
end



parsed_args = parse_args(s)

input_file    = parsed_args["input"]
metadata_file = parsed_args["metadata"]
sensitive     = parsed_args["sensitive"]
heterogeneous = parsed_args["heterogeneous"]
max_k         = parsed_args["max_k"]
n_obs_min     = parsed_args["n_obs_min"]
alpha         = parsed_args["alpha"]


println("hello friend")

# Run FlashWeave with or without metadata
if metadata_file != ""
    net = learn_network(
        input_file, 
        metadata_file,
        sensitive=sensitive,
        heterogeneous=heterogeneous,
        max_k=max_k,
        n_obs_min=n_obs_min,
        alpha=alpha
    )
else
    net = learn_network(
        input_file,
        sensitive=sensitive,
        heterogeneous=heterogeneous,
        max_k=max_k,
        n_obs_min=n_obs_min,
        alpha=alpha
    )
end

save_network("flashweave.edgelist", net)
