#!/usr/bin/env python
import os
import glob
import unittest
import inmembrane


description = """
Inmembrane unit tests. Inmembrane is a proteome annotation pipeline. 
Tests are stored in the `tests` directory and can be run individually.
(c) 2011 Bosco Ho and Andrew Perry

Some of the tests query websites, which may take several minutes:
"""

print description

module_dir = os.path.abspath(os.path.dirname(__file__))
file_tag = os.path.join(module_dir, 'tests', 'test*.py')
test_names = [ \
    os.path.basename(f)[:-3] for f in glob.glob(file_tag)]

for test_name in test_names:
  exec('from tests.%s import *' % test_name)

unittest.main()
