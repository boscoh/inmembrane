import unittest
from semantic_version import Version as SemanticVersion

from inmembrane.tests.PluginTestBase import PluginTestBase
import inmembrane


class TestBomp(PluginTestBase):
    _plugin_name = "bomp_web"

    def test_version(self):
        # parse the version number to ensure it's valid (invalid version number
        # formats will raise a ValueError)
        version = SemanticVersion(inmembrane.__version__)
        self.assertEqual(str(version), inmembrane.__version__)


if __name__ == '__main__':
    unittest.main()
