"""
The network-level steps of the on-the-fly version, run with the network workflow of the pipeline
(pipeline/workflows/network.nf): network inference with FlashWeave, FAPROTAX and clustering with manta.

The server image has the tools and Nextflow (Dockerfile.v2), so the workflow runs with the `local` profile,
without a container engine. MICROBETAG_PIPELINE points to the pipeline; without it, these steps are not available.
"""

import os
import subprocess

import yaml

from microbetag import mtg_logger

logger = mtg_logger(__file__)

PIPELINE = os.environ.get("MICROBETAG_PIPELINE")


def _value(v):
    if isinstance(v, str):
        if v.lower() in ("true", "false"):
            return v.lower() == "true"
        try:
            return int(v)
        except ValueError:
            return v
    return v


def pipeline_available() -> bool:
    return bool(PIPELINE) and os.path.exists(os.path.join(PIPELINE, "main.nf"))


def run_network_workflow(
    out_dir: str, abundance_table: str, network: str = None, metadata: str = None,
    faprotax: bool = False, net_cluster: bool = False, flashweave_args: dict = None,
) -> dict:
    """
    Run the network workflow in out_dir/pipeline.

    Args:
        abundance_table: tab-separated, with the sequence ids in the first column and the taxonomy in the last one
            ('taxonomy'), as the server writes it
        network: a 3-column edge list; None to infer it with FlashWeave

    Returns:
        The outputs: 'edge_list' (the network given or inferred), 'faprotax_sub_tables' and 'manta_net'
        when computed.
    """
    run_dir = os.path.join(out_dir, "pipeline")
    os.makedirs(run_dir, exist_ok=True)
    results = os.path.join(run_dir, "results")

    params = {
        "abundance_file": abundance_table,
        "taxonomy_colname": "taxonomy",
        "delimiter": "\\t",
        "metadata_file": metadata,
        "network": network,
        "faprotax": faprotax,
        "net_cluster": net_cluster,
        "outdir": results,
        # MGG sends the parameters as strings
        "flashweave_args": {k: _value(v) for k, v in (flashweave_args or {}).items()},
    }
    # FlashWeave needs 20 samples by default (n_obs_min); allow fewer, with less reliable tests
    if network is None and "n_obs_min" not in params["flashweave_args"]:
        with open(abundance_table) as f:
            n_samples = len(f.readline().rstrip("\n").split("\t")) - 2
        if n_samples < 20:
            logger.warning(f"Only {n_samples} samples: FlashWeave's tests will be less reliable.")
            params["flashweave_args"]["n_obs_min"] = n_samples

    params_file = os.path.join(run_dir, "params.yaml")
    with open(params_file, "w") as f:
        yaml.safe_dump(params, f, sort_keys=False)

    cmd = [
        "nextflow", "-log", os.path.join(run_dir, "nextflow.log"),
        "run", os.path.join(PIPELINE, "main.nf"), "-profile", "local",
        "-params-file", params_file, "--workflow", "network", "-work-dir", os.path.join(run_dir, "work"),
    ]
    logger.info("Running the network workflow: " + " ".join(cmd))
    proc = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True)
    if proc.returncode != 0:
        logger.error(proc.stdout[-3000:] + proc.stderr[-3000:])
        raise RuntimeError("The network workflow of microbetag-nf failed; see " + run_dir)

    outputs = {"edge_list": network or os.path.join(results, "flashweave", "fw_net.txt")}
    if faprotax:
        outputs["faprotax_sub_tables"] = os.path.join(results, "faprotax", "sub_tables")
    if net_cluster:
        outputs["manta_net"] = os.path.join(results, "manta", "manta_net.cyjs")
    return outputs
