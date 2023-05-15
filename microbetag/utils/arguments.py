"""
Parse arguments
"""

import argparse
from argparse import RawTextHelpFormatter 
import textwrap
import os 

# main_path  = '/'.join(os.getcwd().split("/")[:-1])
# default_config = main_path + "/config.yml"
parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent('''\
                microbetag: annotating microbial co-occurrence networks
                -----------------------------------------------

                microbetag aims to enhance microbial co-occurrence network analysis using 
                data integration and metabolic modeling oriented techniques.
                In case you already have built a network, please provide both the network 
                and its corresponing OTU table files.
                If a network is not availaible, then microbetag will built one for you using 
                FlashWeave and you need to provide is your OTU table, along with its 
                corresponding metadata.
                You should use a configuration file (e.g., config.yml) to pass all your arguments 
                and input files. 
                Further, hints and tips can be found in the template configuration file (config.yml).
                ''')
            )


parser.add_argument('-c', '--configuration_file', dest = 'conf', metavar = '',
                    required = True,
                    help = 'YAML configuration file with the parameter values')

"""
Parse arguments
"""
args = parser.parse_args()
