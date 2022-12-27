"""
Parse arguments
"""

import argparse
import os 

main_path  = '/'.join(os.getcwd().split("/")[:-1])
default_config = main_path + "/config.yml"
parser = argparse.ArgumentParser()

parser.add_argument('-conf', '--configuration_file', metavar ='configuration_file', dest = 'conf',
                    required = True,
                    default = default_config,
                    help = 'YAML configuration file with the parameter values')

"""
Parse arguments
"""
args = parser.parse_args()
