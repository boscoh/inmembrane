import sys, os, math, glob, re, subprocess

default_parms_str = """{
  'organism': 'gram+',
  'signalp4_bin': 'signalp',
  'lipop1_bin': 'LipoP',
  'tmhmm_bin': 'tmhmm',
  'hmmsearch3_bin': 'hmmsearch',
  'memsat3_bin': 'runmemsat',
  'hmm_profiles_dir': '%(hmm_profiles)s',
  'hmm_evalue_cutoff': 0.01,
  'terminal_exposed_loop_min': 50,
  'internal_exposed_loop_min': 100,
}
"""
def get_parms():
  module_dir = os.path.abspath(os.path.dirname(__file__))
  config = os.path.join(module_dir, 'inmembrane.config')
  if not os.path.isfile(config):
    print "# Generating default config", os.path.abspath(config)
    abs_hmm_profiles = os.path.join(module_dir, 'hmm_profiles')
    default_str = default_parms_str % \
        { 'hmm_profiles': abs_hmm_profiles }
    open('inmembrane.config', 'w').write(default_str)
  else:
    print "# Loading existing inmembrane.config"
  parms = eval(open(config).read())
  return parms


def basename(parms):
  return '.'.join(os.path.splitext(parms['fasta'])[:-1])


def run_with_output(cmd):
  p = subprocess.Popen(
      cmd, shell=True, stdout=subprocess.PIPE, 
      stderr=subprocess.PIPE)
  return p.stdout.read()


def run(cmd, out_file="log.txt"):
  if not out_file:
    out_file = "/dev/null"
  binary = cmd.split()[0]
  is_binary_there = False
  if os.path.isfile(binary):
    is_binary_there = True
  if run_with_output('which ' + binary):
    is_binary_there = True
  if not is_binary_there:
    print "# Error: couldn't run executable " + binary
    sys.exit(1)
  full_cmd = cmd + " > " + out_file
  print "#", full_cmd
  if os.path.isfile(out_file) and (out_file != None):
    print "# -> skipped: %s already exists" % out_file
  else:
    os.system(full_cmd)


def read_fasta_keys(fname):
  return [l[1:] for l in open(fname) if l.startswith(">")]


def get_fasta_seq_by_id(fname, prot_id):
  f = open(fname)
  l = f.readline()
  while l:
    if l.startswith(">") and (l[1:].split()[0] == prot_id):
      seq = ""
      l = f.readline()
      while l and not l.startswith(">"):
        seq += l.strip()
        l = f.readline()
      f.close()
      return seq

    l = f.readline()
  f.close()


def hmmsearch3(parms, proteins):
  base = basename(parms)
  profiles_dir = parms['hmm_profiles_dir']

  for hmm_profile in glob.glob(profiles_dir + '/*.hmm'):
    parms['hmm_profile'] = hmm_profile
    hmm_profile = os.path.basename(parms['hmm_profile'])
    hmm_name = hmm_profile.replace('.hmm', '')
    hmmsearch3_out = '%s.hmm.%s.out' % (base, hmm_name)
    run('%(hmmsearch3_bin)s -Z 2000 -E 10 %(hmm_profile)s %(fasta)s' % \
          parms, hmmsearch3_out)
    name = None
    for l in open(hmmsearch3_out):
      words = l.split()
      if l.startswith(">>"):
        name = words[1]
        if 'hmmsearch' not in proteins[name]:
          proteins[name]['hmmsearch'] = []
        continue
      if name is None:
        continue
      if 'conditional E-value' in l:
        evalue = float(words[-1])
        if evalue < parms['hmm_evalue_cutoff']:
          proteins[name]['hmmsearch'].append(hmm_name)


def signalp4(parms, proteins):
  signalp4_out = basename(parms) + '.signalp.out'
  run('%(signalp4_bin)s -t %(organism)s  %(fasta)s' % parms, signalp4_out)
  for l in open(signalp4_out):
    if l.startswith("#"):
      continue
    words = l.split()
    name = words[0]
    proteins[name].update({ 
      'is_signalp': (words[9] == "Y"),
      'signalp_cleave_position': int(words[4]),
    })


def lipop1(parms, proteins):
  lipop1_out = basename(parms) + '.lipop.out' 
  run('%(lipop1_bin)s %(fasta)s' % parms, lipop1_out)
  for l in open(lipop1_out):
    words = l.split()
    if 'score' in l:
      name = words[1]
      if 'cleavage' in l:
        pair = words[5].split("=")[1]
        i = int(pair.split('-')[0])
      else:
        i = None
      proteins[name].update({
        'is_lipop': 'Sp' in words[2],
        'lipop_cleave_position': i,
      })


def tmhmm(parms, proteins):
  tmhmm_out = basename(parms) + '.tmhmm.out'
  run('%(tmhmm_bin)s %(fasta)s' % parms, tmhmm_out)
  name = None
  for l in open(tmhmm_out):
    words = l.split()
    if not words:
      continue
    if l.startswith("#"):
      name = words[1]
    else:
      name = words[0]
    if name is None:
      continue
    if 'tmhmm_helices' not in proteins[name]:
      proteins[name].update({
        'n_tmhmm_helix':0, 
        'sequence_length':0,
        'tmhmm_helices':[],
        'tmhmm_inner_loops':[],
        'tmhmm_outer_loops':[]
      })
    if 'Number of predicted TMHs' in l:
      proteins[name]['n_tmhmm_helix'] = int(words[-1])
    if 'Length' in l:
      proteins[name]['sequence_length'] = int(words[-1])
    if 'inside' in l:
      proteins[name]['tmhmm_inner_loops'].append(
          (int(words[-2]), int(words[-1])))
    if 'outside' in l:
      proteins[name]['tmhmm_outer_loops'].append(
          (int(words[-2]), int(words[-1])))
    if 'TMhelix' in l:
      proteins[name]['tmhmm_helices'].append(
          (int(words[-2]), int(words[-1])))


def tmp_fasta_file(prot_id, seq):
  # FIXME: use proper python temp file, or make a directory
  #        full of output files named by sequence ID
  tmpseq = open('tmp.fasta','w')
  tmpseq.write(">%s\n%s\n" % (prot_id, seq))
  tmpseq.close()
  return 'tmp.fasta'
  
  
def memsat3(parms, proteins):
  """
  Runs MEMSAT3 and parses the output files. Takes a standard 'inmembrane'
  params dictionary and a global proteins dictionary which it populates with
  results.
  
  In the current implementation, this function extracts and feeds sequences to MEMSAT3
  one by one via a temporary file.
  
  These keys are added to the proteins dictionary: 
    'n_memsat3_helix', the number of predicted transmembrane helices;
  
    'memsat3_helices', a list of tuples describing the first and last residue
     number of each transmembrane segment; 
  
    'memsat3_scores', a list of confidence scores (floats) for each predicted 
     tm segment;
  
    'memsat3_inner_loops', a list of tuples describing the first and last residue
     number of each predicted internal loop segment;
  
    'memsat3_outer_loops', a list of tuples describing the first and last residue
     number of each predicted outer loop segment;
  """

  memsat3_out = basename(parms) + '.memsat3.out'
  for prot_id in proteins:
    seq = get_fasta_seq_by_id(parms['fasta'], prot_id)
    sequence_length = len(seq)
    temp_fasta = tmp_fasta_file(prot_id, seq)
    parms['temp_fasta'] = temp_fasta
    run('%(memsat3_bin)s %(temp_fasta)s' % parms, 'tmp.memsat')

    protein = proteins[prot_id]
    protein.update({
      'n_memsat3_helix':0, 
      'sequence_length':len(seq),
      'memsat3_helices':[],
      'memsat3_inner_loops':[],
      'memsat3_outer_loops':[]
    })
    memsat3_helices = protein['memsat3_helices']
    memsat3_scores = protein['memsat3_scores']
    inner_loops = protein['memsat3_inner_loops']
    outer_loops = protein['memsat3_outer_loops']

    f = open('tmp.memsat')
    l = f.readline()
    # parse the final prediction scores from MEMSAT3 output
    while l:
      l = f.readline()
      if l == "FINAL PREDICTION\n":
        num_tms = 0
        f.readline()
        l = f.readline()
        s = l.split(":")
        # parse tm spanning residues and confidence scores
        while re.match("\d", l[0]):
          tokens = s[1].strip().split()
          tok_offset = 0
          if len(tokens) > 2:
            tok_offset = 1
            side_of_membrane_nterminus = tokens[0][1:-1] # 'in' or 'out'
          i = int(tokens[0+tok_offset].split('-')[0])
          j = int(tokens[0+tok_offset].split('-')[1])
          memsat3_helices.append((i, j))
          memsat3_score = float(tokens[1+tok_offset][1:-1])
          memsat3_scores.append(memsat3_score)
          l = f.readline()
          s = l.split(":")
          protein['n_memsat3_helix'] += 1
        f.readline()
        tm_pred = ""
        seq = ""

        # record inner and outer loops
        # TODO: optionally parse out inner and outer helix caps ?
        #       eg OOOO and IIII - are these likely to be cleaved
        #          regions in a shaving experiment ?
        if side_of_membrane_nterminus == 'out':
          current_side = 'out'
        elif side_of_membrane_nterminus == 'in':
          current_side = 'in'
        loop_start = 1
        for tm in memsat3_helices:
          loop_end = tm[0] - 1
          loop = (loop_start, loop_end)
          if current_side == 'out':
            outer_loops.append((loop_start, loop_end))
            current_side = 'in'
          elif current_side == 'in':
            inner_loops.append((loop_start, loop_end))
            current_side = 'out'
          loop_start = tm[1] + 1     
        # capture C-terminal loop segment
        loop_start = tm[1]+1
        loop_end = sequence_length
        if current_side == 'out':
          outer_loops.append((loop_start, loop_end))
        elif current_side == 'in':
          inner_loops.append((loop_start, loop_end))

        """
        # CURRENTLY UNUSED
        # parse tm predictions sequence string
        # eg. +++++OOOOXXXXXXIIIII----- etc
        while l != "":
          tm_pred = tm_pred + f.readline()[:-1]
          seq = seq + f.readline()[:-1]
          l = f.readline()
          l = f.readline()
        """
    f.close()
    if len(memsat3_helices) > 0:
      # print proteins[prot_id]
      # print (seq, tm_pred, memsat3_helices, 
      #     protein['n_memsat3_helix'])
      pass
    else:
      # print (None, None, None, None)
      pass

def chop_nterminal_peptide(protein, i_cut):
  protein['sequence_length'] -= i_cut
  for loop_type in ['tmhmm_outer_loops', 'tmhmm_inner_loops', 'tmhmm_helices']:
    loops = protein[loop_type]
    for i in range(len(loops)):
      j, k = loops[i]
      loops[i] = (j - i_cut, k - i_cut)
  for loop_type in ['tmhmm_outer_loops', 'tmhmm_inner_loops', 'tmhmm_helices']:
    loops = protein[loop_type]
    for i in reversed(range(len(loops))):
      j, k = loops[i]
      # tests if this loop has been cut out
      if j<=0 and k<=0:
        del loops[i]
      # otherewise, neg value means loop is at the new N-terminal
      elif j<=0 and k>0:
        loops[i] = (1, k)
  protein['n_tmhmm_helix'] = len(protein['tmhmm_helices'])


def has_surface_exposed_loop(parms, protein, program='tmhmm'):
  terminal_exposed_loop_min = parms['terminal_exposed_loop_min']
  internal_exposed_loop_min = parms['internal_exposed_loop_min']

  tmhmm_outer_loops = protein['%s_outer_loops' % (program)]
  sequence_length = protein['sequence_length']
  has_no_transmembrane_helices = protein['n_%s_helix' % (program)] == 0

  loop_len = lambda loop: abs(int(loop[1])-int(loop[0])) + 1

  if has_no_transmembrane_helices:
    # treat protein as one entire exposed loop
    return sequence_length >= terminal_exposed_loop_min

  if not tmhmm_outer_loops:
    return False

  # if the N-terminal loop sticks outside
  if tmhmm_outer_loops[0][0] == 1:
    nterminal_loop = tmhmm_outer_loops[0]
    del tmhmm_outer_loops[0]
    if loop_len(nterminal_loop) >= terminal_exposed_loop_min:
      return True

  # if the C-terminal loop sticks outside
  if tmhmm_outer_loops:
    if tmhmm_outer_loops[-1][-1] == sequence_length:
      cterminal_loop = tmhmm_outer_loops[-1]
      del tmhmm_outer_loops[-1]
      if loop_len(cterminal_loop) >= terminal_exposed_loop_min:
        return True

  # test remaining outer loops for length
  for loop in tmhmm_outer_loops:
    if loop_len(loop) >= internal_exposed_loop_min:
      return True

  return False


def print_protein(protein):
  for k in protein.keys():
    # if 'tmh' in k:
      print " %s=%-5s" % (k, protein[k]) 


def predict_surface_exposure(parms, protein):
  if len(protein['hmmsearch']) > 0:
    return "hmmsearch;", "PSE"

  s = ""
  
  if protein['is_lipop']: 
    s += "lipop;"
    chop_nterminal_peptide(protein, protein['lipop_cleave_position'])
    if protein['n_tmhmm_helix'] == 0:
      if protein['sequence_length'] < parms['terminal_exposed_loop_min']:
        return s, "MEMBRANE"
      else:
        return s, "PSE"
  elif protein['is_signalp']:
    s += "signalp;"
    chop_nterminal_peptide(protein, protein['signalp_cleave_position'])
    if protein['n_tmhmm_helix'] == 0:
      return s, "SECRETED"

  if protein['n_tmhmm_helix'] > 0:
    s += "tmhmm;"
    if has_surface_exposed_loop(parms, protein, program='tmhmm'):
      return s, "PSE"
    else:
      return s, "MEMBRANE"

  if protein['n_memsat3_helix'] > 0:
    s += "memsat3;"
    if has_surface_exposed_loop(parms, protein, program='memsat3'):
      return s, "PSE"
    else:
      return s, "MEMBRANE"

  return s, "CYTOPLASM"


def identify_pse_proteins(parms, fasta):
  parms['fasta'] = fasta
  # initialize the proteins data structure
  headers = read_fasta_keys(fasta)
  proteins = {}
  prot_ids = []
  for header in headers:
    prot_id = header.split()[0]
    prot_ids.append(prot_id)
    proteins[prot_id] = {
      'hmmsearch': [],
      'is_signalp': False,
      'is_lipop': None,
      'n_tmhmm_helix': 0,
      'name': ' '.join(header.split()[1:]),
    }
  for extract_protein_feature in [memsat3, signalp4]:\
      # FIXME: temporarily commented for testing
      #[signalp4, lipop1, tmhmm, hmmsearch3]:
    extract_protein_feature(parms, proteins)
  for prot_id in prot_ids:
    details, category = \
        predict_surface_exposure(parms, proteins[prot_id])
    if details.endswith(';'):
      details = details[:-1]
    if details is '':
      details = "."
    proteins[prot_id]['details'] = details
    proteins[prot_id]['category'] = category
        
  return prot_ids, proteins



if __name__ == "__main__":
  fasta = sys.argv[1]
  parms = get_parms()
  prot_ids, proteins = identify_pse_proteins(
      parms, fasta)
  for prot_id in prot_ids:
    protein = proteins[prot_id]
    print "%-15s %-13s %-20s %s" % \
        (prot_id.split("|")[1], 
         protein['category'], 
         protein['details'],
         protein['name'][:60])


  





