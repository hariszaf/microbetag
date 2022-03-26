"""
Parse arguments
"""

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-conf', '--configuration_file', metavar ='configuration_file', dest = 'conf',
                    required = True,
                    help='YAML configuration file with the parameter values')

"""
Parse arguments
"""
args = parser.parse_args()
