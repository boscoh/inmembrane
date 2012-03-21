#
# Common helper functions for inmembrane
#

import os, subprocess, sys



LOG_DEBUG = False
LOG_SILENT = False

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
  log_stderr("# " + cmd + " > " + out_file)
  if os.path.isfile(out_file) and (out_file != None):
    log_stderr("# -> skipped: %s already exists" % out_file)
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


def silence_log(b):
  global LOG_SILENT
  LOG_SILENT = b


def log_stderr(s):
  if LOG_SILENT:
    return
  if s and s[-1] != "\n":
    s += "\n"
  if not s.startswith("#"):
    s = "#" + s
  sys.stderr.write(s)


def log_stdout(s):
  if LOG_SILENT:
    return
  print s


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


def create_proteins_dict(fasta):
  """
  The main data-structure used in inmembrane is the proteins dictionary,
  which is initialized here from a source FASTA file. Also returns a list 
  of seqID's in the same order as that in the source FASTA file.
  """
  seqids = []
  seqid = None
  proteins = {}
  for l in open(fasta):
    if l.startswith(">"):
      seqid, name = parse_fasta_header(l)
      seqids.append(seqid)
      proteins[seqid] = {
        'seq':"",
        'name':name,
      }
      continue
    if seqid is not None:
      words = l.split()
      if words:
        proteins[seqid]['seq'] += words[0]
  return seqids, proteins
  

def chop_nterminal_peptide(protein, i_cut):
  protein['sequence_length'] -= i_cut
  for prop in protein:
    if '_loops' in prop or '_helices' in prop:
      loops = protein[prop]
      for i in range(len(loops)):
        j, k = loops[i]
        loops[i] = (j - i_cut, k - i_cut)
  for prop in protein:
    if '_loops' in prop or '_helices' in prop:
      loops = protein[prop]
      for i in reversed(range(len(loops))):
        j, k = loops[i]
        # tests if this loop has been cut out
        if j<=0 and k<=0:
          del loops[i]
        # otherewise, neg value means loop is at the new N-terminal
        elif j<=0 and k>0:
          loops[i] = (1, k)


def print_proteins(proteins):
  print "{"
  for seqid in proteins:
    print "  '%s': {" % seqid
    for key, value in proteins[seqid].items():
      print "    '%s': %s, " % (key, repr(value))
    print "  },"
  print "}"

  
def eval_surface_exposed_loop(
    sequence_length, n_transmembrane_region, outer_loops, 
    terminal_exposed_loop_min, internal_exposed_loop_min):
    
  if n_transmembrane_region == 0:
    # treat protein as one entire exposed loop
    return sequence_length >= terminal_exposed_loop_min

  if not outer_loops:
    return False

  loop_len = lambda loop: abs(loop[1]-loop[0]) + 1

  # if the N-terminal loop sticks outside
  if outer_loops[0][0] == 1:
    nterminal_loop = outer_loops[0]
    del outer_loops[0]
    if loop_len(nterminal_loop) >= terminal_exposed_loop_min:
      return True

  # if the C-terminal loop sticks outside
  if outer_loops:
    if outer_loops[-1][-1] == sequence_length:
      cterminal_loop = outer_loops[-1]
      del outer_loops[-1]
      if loop_len(cterminal_loop) >= terminal_exposed_loop_min:
        return True

  # test remaining outer loops for length
  for loop in outer_loops:
    if loop_len(loop) >= internal_exposed_loop_min:
      return True

  return False


def clean_directory(top, excluded_files):
  """
  Deletes an entire directory tree, equivalent to rm -rf <top>, 
  with an excluded file list.
  """
  for root, dirs, files in os.walk(top, topdown=False):
    for name in files:
      if name not in excluded_files:
        os.remove(os.path.join(root, name))
    for name in dirs:
      os.rmdir(os.path.join(root, name))


