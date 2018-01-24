import os
import glob
import unittest
import sys

import inmembrane
import inmembrane.tests

from inmembrane import helpers
from inmembrane.plugins import bomp_web

from inmembrane.tests.PluginTestBase import PluginTestBase


class TestBomp(PluginTestBase):
    _plugin_name = "bomp_web"

    def test_bomp(self):
        self.expected_output = {
            u'gi|107837101': 3,
            u'gi|107836588': 5,
            u'gi|107836852': 5
        }

        self.output = bomp_web.annotate(self.params, self.proteins, force=True)

        self.assertEqual(self.expected_output, self.output)


if __name__ == '__main__':
    unittest.main()
