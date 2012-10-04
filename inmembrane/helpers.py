#
# Common helper functions for inmembrane
#

import inmembrane
import os, subprocess, sys
import textwrap

LOG_SILENT = False

def dict_get(this_dict, prop):
  """
  Useful helper function to get values from dicts. Takes in
  the possibility that the key may not exist, and in that
  case returns False. Makes it easier to write code to avoid
  handling the case of missing keys.
  """
  if prop not in this_dict:
    return False
  return this_dict[prop]
  

def run_with_output(cmd):
  """
  Runs an external program as a child process and captures
  the input. May not work with external programs that include
  complicated levels of shells.
  """
  p = subprocess.Popen(
      cmd, shell=True, stdout=subprocess.PIPE, 
      stderr=subprocess.PIPE)
  return p.stdout.read()


def run(cmd, out_file=None):
  """
  Wrapper function to run external program so that existence
  of binary can be checked and the output directed to a specific 
  file.
  """
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
  """
  Turns on logging silent mode w.r.t. to log_stderr and log_stdout
  """
  global LOG_SILENT
  LOG_SILENT = b


def log_stderr(s):
  """
  Wrapper for all stderr out. Allows future customization.
  """
  if LOG_SILENT:
    return
  if s and s[-1] != "\n":
    s += "\n"
  if not s.startswith("#"):
    s = "# " + s
  sys.stderr.write(s)


def log_stdout(s):
  """
  Wrapper for all stderr out. Allows future customization.
  """
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
    name = tokens[-1:][0].strip()
  # otherwise just split on spaces & hope for the best
  else:
    tokens = header.split()
    seqid = tokens[0]
    name = header[0:-1].strip()
  
  return seqid, name
  

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
    proteins[seqid]['sequence_length'] = len(proteins[seqid]['seq'])
  return seqids, proteins
  

def print_proteins(proteins):
  """
  Prints out the proteins dictionary in a formatted 
  manner that is also Python-eval compatible.
  """
  print "{"
  for seqid in proteins:
    print "  '%s': {" % seqid
    for key, value in proteins[seqid].items():
      print "    '%s': %s, " % (key, repr(value))
    print "  },"
  print "}"
  # Standard Library alternative 
  # import pprint
  # pp = pprint.PrettyPrinter(indent=4)
  # print pp.pformat(proteins)
 
  
def write_proteins_fasta(
    fasta_filename, proteins, seqids, width=50):
  """
  Creates a fasta file of the sequences of a subset of the proteins.
  """
  f = open(fasta_filename, "w")
  for seqid in seqids:
    seq_wrap = textwrap.fill(proteins[seqid]['seq'], width)
    f.write(">%s\n%s\n" % (proteins[seqid]['name'], seq_wrap))
  f.close()


def chop_nterminal_peptide(protein, i_cut):
  """
  Finds all transmembrane fields in protein ("*_loops" and "*_helices")
  and deletes an N-terminal section of the protein indicated by i_cut.

  Assuming the topology of inner and outer loops remain the same, this
  may delete a certain number of elements at the N-terminus. Allows a
  simple way of removing lipo-protein and secretion-protein signals from
  the evaluation of the outer_loops of the protein.
  """
  protein['sequence_length'] -= i_cut
  for prop in protein:
    if '_loops' in prop or '_helices' in prop:
      sses = protein[prop]
      for i in range(len(sses)):
        j, k = sses[i]
        sses[i] = (j - i_cut, k - i_cut)
  for prop in protein:
    if '_loops' in prop or '_helices' in prop:
      sses = protein[prop]
      for i in reversed(range(len(sses))):
        j, k = sses[i]
        # tests if this loop or TM-helix has been cut out
        if j<=0 and k<=0:
          del sses[i]
        # otherewise, neg value means loop or TM-helix is at the new N-terminal
        elif j<=0 and k>0:
          sses[i] = (1, k)
          # if a TM-helix gets cut in half and becomes a new N-terminal,
          # convert the remaining residues to a loop and remove the helix
          if '_helices' in prop:
            program = prop.split('_')[0]
            for x in protein:
              if x == '%s_loops' % program:
                new_N_loop = protein[x][0]
                new_N_loop[0] = 1
            del sses[i]
            

def eval_surface_exposed_loop(
    sequence_length, n_transmembrane_region, outer_loops, 
    terminal_exposed_loop_min, internal_exposed_loop_min):
  """
  This is the key algorithm in SurfG+ to identify Potentially Surface
  Exposed proteins. It evaluates all loops that poke out of the periplasmic
  side of the Gram+ bacterial membrane and tests if the loops are long
  enough to stick out of the peptidoglycan layer to be cleaved by proteases
  in a cell-shaving experiment.

  Returns True if any outer loop, or the N- or C-terminii are
  longer than the given thresholds.
  """

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
