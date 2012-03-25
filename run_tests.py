#!/usr/bin/env python
import sys, os
import glob
import unittest
import inmembrane
import helpers

description = """
Inmembrane unit tests. Inmembrane is a proteome annotation pipeline. 
Tests are stored in the `tests` directory and can be run individually.
(c) 2011 Bosco Ho and Andrew Perry

Some of the tests query websites, which may take several minutes:
"""

print description

from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-n", "--no-network",
                     action="store_true", dest="no_network", default=False,
                     help="don't run tests that require a network connection")
(options, args) = optparser.parse_args()
# unusual hack to remove our custom flags, since unittest expects it's own set
if "-n" in sys.argv:
  sys.argv.remove("-n")
if "--no-network" in sys.argv:
  sys.argv.remove("--no-network")
  
module_dir = os.path.abspath(os.path.dirname(__file__))
file_tag = os.path.join(module_dir, 'tests', 'test*.py')
test_names = [ \
    os.path.basename(f)[:-3] for f in glob.glob(file_tag)]

tests_to_run = []
for test_name in test_names:
  if not (options.no_network and (test_name[-4:] == "_web")):
    exec('from tests.%s import *' % test_name)
    tests_to_run.append(test_name)

helpers.log_stderr("Will run tests for: " + ",".join(tests_to_run))

unittest.main()
