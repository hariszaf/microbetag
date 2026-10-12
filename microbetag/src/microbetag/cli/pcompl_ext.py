"""Extend pathway complementarities with module descriptions (microbetag-pcompl-ext)."""
import argparse
import os
from pathlib import Path

from microbetag.cli import KEGG_MAPPINGS


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-pcompl-ext", description=__doc__)
    parser.add_argument("pc_file")
    parser.add_argument("pc_percent", type=float)
    parser.add_argument("threads", type=int)
    args = parser.parse_args(argv)

    from microbetag.utils import extend_complements

    extend_complements(
        complements_json = args.pc_file,
        descrps_path     = os.path.join(KEGG_MAPPINGS, "module_descriptions"),
        pc_percent       = args.pc_percent,
        pc_dir           = Path(args.pc_file).parent.absolute(),
        n_workers        = args.threads,
    )


if __name__ == "__main__":
    main()
