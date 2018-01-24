import os, tempfile, shutil
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers


class PluginTestBase(unittest.TestCase):
    _plugin_name = ""

    def setUp(self):
        """
        Sets up a directory, a parameters dictionary and a proteins dictionary
        object in preparation for running a test of a plugin.

        Creates a temporary directory (eg /tmp/.inmembrane_signalp_web_TleeRw ),
        copies the test input data (input.fasta) into it and changes the
        current working directory.

        Should be subclassed adding a method with the specific test(s) for
        the plugin. The subclass should set self._plugin_name to the name of
        the plugin.
        """
        self.test_data_dir = os.path.join(
            os.path.abspath(
                os.path.dirname(inmembrane.tests.__file__)), self._plugin_name)

        self.output_dir = tempfile.mkdtemp(
            prefix=".inmembrane_%s_" % (self._plugin_name))
        shutil.copyfile(os.path.join(self.test_data_dir, "input.fasta"),
                        os.path.join(self.output_dir, "input.fasta"))
        os.chdir(self.output_dir)
        helpers.silence_log(True)

        self.params = inmembrane.get_params()
        self.params['fasta'] = "input.fasta"
        self.seqids, self.proteins = \
            helpers.create_proteins_dict(self.params['fasta'])
