import getopt
import subprocess
import os, sys, time



opts, args = getopt.getopt(sys.argv[1:], "w:", ["workflow"])

# Run the create-yaml.py script
subprocess.call(['python3','scripts/create-yaml.py', args[0]])


# Once the .yml is ready, run the microbetag workflow
subprocess.call(["cwl-runner", "--debug", "microbetag.cwl", "microbetag-job.yml"])



