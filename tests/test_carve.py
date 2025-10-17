import os
import yaml
import unittest
from pathlib import Path

from microbetag.config import Config
#from microbetag.genres import GEMSReconstruction
from microbetag.wrappers import build_genres

root_dir  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
test_data = os.path.join(root_dir, "test_data", "test_carve")

# Get configuration based on the config YAML file
config_file = os.path.join(test_data, "config_carve.yml")
with open(config_file, "r") as yaml_file:
    config = Config(yaml.safe_load(yaml_file), config_file)

# Remove any .xml file that may be in the output folder from previous tests
for file_path in Path(config.genres).glob("*.xml"):
    file_path.unlink()  # This deletes the file


class testBuildGemWithCarve(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        cls.config = config
        cls.GEMSReconstruction = None

    def test_carve(self):

        build_genres(self.config)


if __name__ == "__main__":
    unittest.main()
