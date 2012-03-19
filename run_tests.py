#!/usr/bin/env python
import os
import glob
import unittest
import inmembrane


module_dir = os.path.abspath(os.path.dirname(__file__))
file_tag = os.path.join(module_dir, 'tests', 'test*.py')
test_names = [ \
    os.path.basename(f)[:-3] for f in glob.glob(file_tag)]

# test_names = ['test_hmmsearch3']

for test_name in test_names:
  exec('from tests.%s import *' % test_name)

unittest.main()
