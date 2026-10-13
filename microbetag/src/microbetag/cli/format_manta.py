"""Convert a 3-column edgelist to a .cyjs manta input (microbetag-format-manta)."""
import argparse
from pathlib import Path


class Config:
    def __init__(self, abundance_file, network):
        self.abundance_file    = abundance_file
        self.network           = network
        # Name of the formatted network file (output)
        self.base_network_file = Path(network).with_suffix(".cyjs")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="microbetag-format-manta", description=__doc__)
    parser.add_argument("abundance_file")
    parser.add_argument("network")
    args = parser.parse_args(argv)

    from microbetag.cli.annotate import load_data
    from microbetag.helpers import manta_input_net

    manta_input_net(load_data(Config(args.abundance_file, args.network)))


if __name__ == "__main__":
    main()
