#
# Common helper functions for inmembrane
#

import os, subprocess


__DEBUG__ = False


def dict_get(this_dict, prop):
  if prop not in this_dict:
    return False
  return this_dict[prop]
  

def run_with_output(cmd):
  p = subprocess.Popen(
      cmd, shell=True, stdout=subprocess.PIPE, 
      stderr=subprocess.PIPE)
  return p.stdout.read()


def run(cmd, out_file=None):
  error_output("# " + cmd + " > " + out_file)
  if os.path.isfile(out_file) and (out_file != None):
    error_output("# -> skipped: %s already exists" % out_file)
    return
  binary = cmd.split()[0]
  is_binary_there = False
  if os.path.isfile(binary):
    is_binary_there = True
  if run_with_output('which ' + binary):
    is_binary_there = True
  if not is_binary_there:
    raise IOError("Couldn't find executable binary '" + binary + "'")
  if out_file:
    os.system(cmd + " > " + out_file)
  else:
    os.system(cmd)


def error_output(s):
  if not __DEBUG__:
    return
  if s and s[-1] != "\n":
    s += "\n"
  sys.stderr.write(s)


def parse_fasta_header(header):
  """
  Parses a FASTA format header (with our without the initial '>') and returns a
  tuple of sequence id and sequence name/description.
  
  If NCBI SeqID format (gi|gi-number|gb|accession etc, is detected
  the first id in the list is used as the canonical id (see see
  http://www.ncbi.nlm.nih.gov/books/NBK21097/#A631 ).
  """
  # check to see if we have an NCBI-style header
  if header[0] == '>':
    header = header[1:]
  if header.find("|") != -1:
    tokens = header.split('|')
    # "gi|ginumber|gb|accession bla bla" becomes "gi|ginumber"
    seqid = "%s|%s" % (tokens[0], tokens[1].split()[0])
    desc = tokens[-1:][0].strip()
  # otherwise just split on spaces & hope for the best
  else:
    tokens = header.split()
    seqid = tokens[0]
    desc = header[0:-1].strip()
  
  return seqid, desc
  

def seqid_to_filename(seqid):
  """
  Makes a sequence id filename friendly.
  (eg, replaces '|' with '_')
  """
  return seqid.replace("|", "_")
  
