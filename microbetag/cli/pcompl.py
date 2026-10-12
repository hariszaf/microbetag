"""Export pathway complementarities from a merged KO file (microbetag-pcompl)."""
import argparse
from pathlib import Path

from microbetag.cli import KEGG_MAPPINGS


class Config:
    kegg_mappings                = Path(KEGG_MAPPINGS)
    kegg_modules_to_maps         = kegg_mappings / "module_map_pairs.tsv"
    ref_ko_per_module            = kegg_mappings / "kegg_terms_per_module.tsv"
    modules_definitions_json_map = kegg_mappings / "module_definition_map.json"
    module_descriptions          = kegg_mappings / "module_descriptions"

    def __init__(self, alts_file, pc_file, ko_merged, tinyurl, threads):
        self.alts_file = alts_file
        self.pc_file   = pc_file
        self.ko_merged = ko_merged
        self.tinyurl   = tinyurl
        self.threads   = threads


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-pcompl", description=__doc__)
    parser.add_argument("alts_file", help="Output: module alternatives JSON")
    parser.add_argument("pc_file", help="Output: pathway complementarities JSON")
    parser.add_argument("ko_merged", help="3-column file with KO terms per genome")
    parser.add_argument("tinyurl")
    parser.add_argument("threads", type=int)
    args = parser.parse_args(argv)

    from microbetag.pathway_complementarity import export_pathway_complementarities
    from microbetag.utils import load_merged_ko_file

    config   = Config(args.alts_file, args.pc_file, args.ko_merged, args.tinyurl, args.threads)
    pivot_df = load_merged_ko_file(config.ko_merged)
    export_pathway_complementarities(config, pivot_df)


if __name__ == "__main__":
    main()
