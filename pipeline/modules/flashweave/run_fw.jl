using ArgParse
using JSON
using FlashWeave

# Parse arguments
settings = ArgParseSettings()

@add_arg_table settings begin

    "--input"
        help = "Input abundance table"
        arg_type = String
        required = true

    "--metadata"
        help = "Optional metadata file"
        arg_type = String
        default = ""


    "--fw_args"
        help = "FlashWeave args as JSON string"
        arg_type = String
        required = false
        default = "{}"
end


default_fw = Dict(
    "sensitive"        =>  true,
    "heterogeneous"    =>  false,
    "transposed"       =>  true,    # Sos!
    "alpha"            =>  0.01,
    "max_k"            =>  3,
    "track_rejections" =>  false,
    "normalize"        =>  true,
    "time_limit"       =>  -1.0,
    "n_obs_min"        =>  -1,
    "FDR"              =>  true,
    "hps"              =>  5,
    "max_tests"        =>  Int(10e6),
    "feed_forward"     =>  true,
    "conv"             =>  0.01,
    "fast_elim"        =>  true,
    "prec"             =>  32,
    "update_interval"  =>  30.0
)

parsed_args = parse_args(settings)
fw_args     = JSON.parse(parsed_args["fw_args"])  # Dict{String,Any}

# override only provided keys
for (k,v) in fw_args
    if v !== nothing
        default_fw[k] = v
    end
end

# optional: convert to NamedTuple for splatting
fw_kwargs = (; (Symbol(k) => v for (k,v) in default_fw)...)

# Abundance table with ASVs, OTUs, bins or MAGs as rows, and samples as columns!!!
input_file = parsed_args["input"]


# Run FlashWeave with or without metadata
if haskey(parsed_args, "metadata") && parsed_args["metadata"] !== nothing && parsed_args["metadata"] != ""
    metadata_file = parsed_args["metadata"]
    net = learn_network(
        input_file, 
        metadata_file;  # (; ...) → creates a NamedTuple suitable for keyword argument splatting -- fw_kwargs... after ; → passes as keywords, not positionals
        fw_kwargs...
    )
else
    net = learn_network(
        input_file;
        fw_kwargs...
    )
end

save_network("flashweave.edgelist", net)
