import sys
import os
import math
import glob
import re
import subprocess
import shutil
from optparse import OptionParser


default_params_str = """{
  'fasta': '',
  'out_dir': '',
  'organism': 'gram+',
  'signalp4_bin': 'signalp',
  'lipop1_bin': 'LipoP',
  'tmhmm_bin': 'tmhmm',
  'helix_programs': ['tmhmm', 'memsat3'],
  'memsat3_bin': 'runmemsat',
  'hmmsearch3_bin': 'hmmsearch',
  'hmm_profiles_dir': '%(hmm_profiles)s',
  'hmm_evalue_max': 0.1,
  'hmm_score_min': 10,
  'terminal_exposed_loop_min': 50,
  'internal_exposed_loop_min': 100,
}
"""


def get_params():
  module_dir = os.path.abspath(os.path.dirname(__file__))
  config = os.path.join(module_dir, 'inmembrane.config')
  if not os.path.isfile(config):
    print "# Couldn't find inmembrane.config file"
    print "# So, will generate a default config", os.path.abspath(config)
    abs_hmm_profiles = os.path.join(module_dir, 'hmm_profiles')
    default_str = default_params_str % \
        { 'hmm_profiles': abs_hmm_profiles }
    open('inmembrane.config', 'w').write(default_str)
  else:
    print "# Loading existing inmembrane.config"
  params = eval(open(config).read())
  return params


def dict_prop_truthy(this_dict, prop):
  if prop not in this_dict:
    return False
  return this_dict[prop]
  

def run_with_output(cmd):
  p = subprocess.Popen(
      cmd, shell=True, stdout=subprocess.PIPE, 
      stderr=subprocess.PIPE)
  return p.stdout.read()


def run(cmd, out_file=None):
  full_cmd = cmd + " > " + out_file
  print "#", full_cmd
  if os.path.isfile(out_file) and (out_file != None):
    print "# -> skipped: %s already exists" % out_file
    return
  if not out_file:
    out_file = "/dev/null"
  binary = cmd.split()[0]
  is_binary_there = False
  if os.path.isfile(binary):
    is_binary_there = True
  if run_with_output('which ' + binary):
    is_binary_there = True
  if not is_binary_there:
    print "# Error: couldn't find executable " + binary
    sys.exit(1)
  os.system(full_cmd)


def create_protein_data_structure(fasta):
  prot_ids = []
  prot_id = None
  proteins = {}
  for l in open(fasta):
    if l.startswith(">"):
      tokens = l.split()
      prot_id = tokens[0][1:]
      prot_ids.append(prot_id)
      proteins[prot_id] = {
        'seq':"",
        'name':l[1:-1],
      }
      continue
    if prot_id is not None:
      words = l.split()
      if words:
        proteins[prot_id]['seq'] += words[0]
  return prot_ids, proteins


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


def hmmsearch3(params, proteins):
  file_tag = os.path.join(params['hmm_profiles_dir'], '*.hmm')
  for hmm_profile in glob.glob(file_tag):
    params['hmm_profile'] = hmm_profile
    hmm_profile = os.path.basename(params['hmm_profile'])
    hmm_name = hmm_profile.replace('.hmm', '')
    hmmsearch3_out = 'hmm.%s.out' % hmm_name
    run('%(hmmsearch3_bin)s -Z 2000 -E 10 %(hmm_profile)s %(fasta)s' % \
          params, hmmsearch3_out)
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
        score = float(words[-5])
        if evalue <= params['hmm_evalue_max'] and \
            score >= params['hmm_score_min']:
          print evalue, score
          proteins[name]['hmmsearch'].append(hmm_name)


def signalp4(params, proteins):
  signalp4_out = 'signalp.out'
  run('%(signalp4_bin)s -t %(organism)s  %(fasta)s' % params, signalp4_out)
  for l in open(signalp4_out):
    if l.startswith("#"):
      continue
    words = l.split()
    name = words[0]
    proteins[name].update({ 
      'is_signalp': (words[9] == "Y"),
      'signalp_cleave_position': int(words[4]),
    })


def lipop1(params, proteins):
  lipop1_out = 'lipop.out'
  run('%(lipop1_bin)s %(fasta)s' % params, lipop1_out)
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

def bomp_web(parms, proteins, \
             url="http://services.cbu.uib.no/tools/bomp/", force=False):
  """
  Uses the BOMP web service (http://services.cbu.uib.no/tools/bomp/) to
  predict if proteins are outer membrane beta-barrels.
  """
  
  bomp_out = basename(parms) + '.bomp.out'
  print "# BOMP(web) %s > %s" % (parms['fasta'], bomp_out)
  
  if not force and os.path.isfile(bomp_out):
    print "# -> skipped: %s already exists" % bomp_out
    return
  
  import time, StringIO
  import twill
  from twill.commands import find, formfile, fv, go, show, showlinks, submit

  # dump extraneous output here so we don't see it
  twill.set_output(StringIO.StringIO())
  
  go(url)
  #showforms()
  formfile("1", "queryfile", parms["fasta"])
  submit()
  
  # extract the job id from the page
  links = showlinks()
  job_id = None
  for l in links:
    if l.url.find("viewOutput") != -1:
      # grab job id from "viewOutput?id=16745338"
      job_id = int(l.url.split("=")[1])
  
  # print "BOMP job id: ", job_id
  
  if not job_id:
    # something went wrong
    sys.stderr.write("# BOMP error: Can't find job id")
    return
  
  # parse the HTML table and extract categories
  go("viewOutput?id=%i" % (job_id))
  
  polltime = 10
  sys.stderr.write("# Waiting for BOMP to finish .")
  while True:
    try:
      find("Not finished")
      sys.stderr.write(".")
    except:
      # Finished ! Pull down the result page.
      sys.stderr.write(". done!\n")
      go("viewOutput?id=%i" % (job_id))
      # print show()
      break
      
    # Not finished. We keep polling for a time until
    # we give up
    time.sleep(polltime)
    polltime = polltime * 2
    if polltime >= 7200: # 2 hours
      sys.stderr.write("# BOMP error: Taking too long.")
      return
    go("viewOutput?id=%i" % (job_id))
    #print show()
      
  bomp_html = show()
  #print bomp_html
  
  # Results are in the only <table> on this page, formatted like:
  # <tr><th>gi|107836852|gb|ABF84721.1<th>5</tr>
  from BeautifulSoup import BeautifulSoup
  soup = BeautifulSoup(bomp_html)
  bomp_categories = {} # dictionary of {name, category} pairs
  for tr in soup.findAll('tr')[1:]:
    n, c = tr.findAll('th')
    name = n.text.split()[0].strip()
    category = int(c.text)
    bomp_categories[name] = category
  
  # write BOMP results to a tab delimited file
  fh = open(bomp_out, 'w')
  for k,v in bomp_categories.iteritems():
    fh.write("%s\t%i\n" % (k,v))
  fh.close()
  
  #print bomp_categories
  
  # label proteins with bomp classification (int) or False
  for name in proteins:
    if "bomp" not in proteins[name]:
      if name in bomp_categories:
        category = bomp_categories[name]
        proteins[name]['bomp'] = category
      else:
        proteins[name]['bomp'] = False
  
  #print proteins
  
  return bomp_categories
  
  """
  # Alternative: just get binary classification results via the
  #              FASTA output BOMP links to
  #
  # use the job id to jump straight to the fasta results
  # if a sequence is here, it's classified as an OMP barrel
  go("viewFasta?id=%i" % (job_id))
  bomp_seqs = show()
  bomp_fasta_headers = read_fasta_keys(StringIO.StringIO(show()))
  # label the predicted TMBs
  for name in bomp_fasta_headers:
    proteins[name]['bomp'] = True
    
  # label all the non-TMBs
  for name in proteins:
    if "bomp" not in proteins[name]:
      proteins[name]['bomp'] = False
  """


def tmhmm(params, proteins):
  tmhmm_out = 'tmhmm.out'
  run('%(tmhmm_bin)s %(fasta)s' % params, tmhmm_out)
  name = None
  for i_line, l in enumerate(open(tmhmm_out)):
    if i_line == 0:
      continue
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
        'sequence_length':0,
        'tmhmm_helices':[],
        'tmhmm_inner_loops':[],
        'tmhmm_outer_loops':[]
      })
    if 'Number of predicted TMHs' in l:
      n_helix = int(words[-1])
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


def has_transmembrane_in_globmem(globmem_out):
  for l in open(globmem_out):
    if "Your protein is probably not a transmembrane protein" in l:
      return False
  return True


def parse_memsat(protein, memsat_out):
    # parse tm spanning residues and confidence scores
    f = open(memsat_out)
    l = f.readline()
    while l:
      l = f.readline()
      if l == "FINAL PREDICTION\n":
        f.readline()
        l = f.readline()
        s = l.split(":")
        while re.match("\d", l[0]):
          tokens = s[1].strip().split()
          tok_offset = 0
          if len(tokens) > 2:
            tok_offset = 1
            side_of_membrane_nterminus = tokens[0][1:-1] # 'in' or 'out'
          i = int(tokens[tok_offset].split('-')[0])
          j = int(tokens[tok_offset].split('-')[1])
          protein['memsat3_helices'].append((i, j))
          score = float(tokens[1+tok_offset][1:-1])
          protein['memsat3_scores'].append(score)
          l = f.readline()
          s = l.split(":")
        f.readline()
        
        # record inner and outer loops
        inner_loops = protein['memsat3_inner_loops']
        outer_loops = protein['memsat3_outer_loops']
        sequence_length = protein['sequence_length']
        if side_of_membrane_nterminus == 'out':
          loops = outer_loops
        elif side_of_membrane_nterminus == 'in':
          loops = inner_loops
        loop_start = 1
        for tm in protein['memsat3_helices']:
          loop_end = tm[0] - 1
          loop = (loop_start, loop_end)
          loops.append((loop_start, loop_end))
          if loops == outer_loops:
            loops = inner_loops
          else:
            loops = outer_loops
          loop_start = tm[1] + 1     
        # capture C-terminal loop segment
        loop_start = tm[1]+1
        loop_end = sequence_length
        loops.append((loop_start, loop_end))

    f.close()

    
    
def memsat3(params, proteins):
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
  for prot_id in proteins:
    protein = proteins[prot_id]
    seq = protein['seq']
    protein.update({
      'sequence_length':len(seq),
      'memsat3_scores':[],
      'memsat3_helices':[],
      'memsat3_inner_loops':[],
      'memsat3_outer_loops':[]
    })

    # write seq to single fasta file
    single_fasta = prot_id.replace("|", "_") + '.fasta'
    if not os.path.isfile(single_fasta):
      open(single_fasta, 'w').write(">%s\n%s\n" % (prot_id, seq))
    memsat_out = single_fasta.replace('fasta', 'memsat')
    run('%s %s' % (params['memsat3_bin'], single_fasta), memsat_out)

    globmem_out = single_fasta.replace('fasta', 'globmem')
    if has_transmembrane_in_globmem(globmem_out):
      parse_memsat(protein, memsat_out)
      

def chop_nterminal_peptide(protein, i_cut):
  protein['sequence_length'] -= i_cut
  print "#"
  print "# cut %s:" % protein['name'].strip()
  for prop in protein:
    if '_loops' in prop or '_helices' in prop:
      loops = protein[prop]
      print "# original %s: %s" % (prop, loops)
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
      print "# cut at %d %s: %s" % (i_cut, prop, loops)


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


def predict_surface_exposure(params, protein):

  def sequence_length(protein):
    return protein['sequence_length']
    
  def has_tm_helix(protein):
    for program in params['helix_programs']:
      if dict_prop_truthy(protein, '%s_helices' % program):
        return True
    return False

  def has_surface_exposed_loop(protein):
    for program in params['helix_programs']:
      if eval_surface_exposed_loop(
          protein['sequence_length'], 
          len(protein['%s_helices' % (program)]), 
          protein['%s_outer_loops' % (program)], 
          params['terminal_exposed_loop_min'], 
          params['internal_exposed_loop_min']):
        return True
    return False

  terminal_exposed_loop_min = \
      params['terminal_exposed_loop_min']

  is_hmm_profile_match = dict_prop_truthy(protein, 'hmmsearch')
  is_lipop = dict_prop_truthy(protein, 'is_lipop')
  if is_lipop:
    i_lipop_cut = protein['lipop_cleave_position']
  is_signalp = dict_prop_truthy(protein, 'is_signalp')
  if is_signalp:
    i_signalp_cut = protein['signalp_cleave_position']

  details = ""
  if is_hmm_profile_match:
    details += "hmm(%s);" % protein['hmmsearch'][0]
  if is_lipop: 
    details += "lipop;"
  if is_signalp:
    details += "signalp;"
  for program in params['helix_programs']:
    if has_tm_helix(protein):
      n = len(protein['%s_helices' % program])
      details += program + "(%d);" % n

  if is_lipop: 
    chop_nterminal_peptide(protein, i_lipop_cut)
  elif is_signalp:
    chop_nterminal_peptide(protein, i_signalp_cut)

  if is_hmm_profile_match:
    category =  "PSE"
  elif has_tm_helix(protein):
    if has_surface_exposed_loop(protein):
      category = "PSE"
    else:
      category = "MEMBRANE"
  else:
    if is_lipop:
      # whole protein considered outer terminal loop
      if sequence_length(protein) < terminal_exposed_loop_min:
        category = "MEMBRANE"
      else:
        category = "PSE"
    elif is_signalp:
      category = "SECRETED"
    else:
      category = "CYTOPLASM"

  return details, category

def identify_pse_proteins(params):
  if dict_prop_truthy(params, 'out_dir'):
    base_dir = params['out_dir']
  else:
    base_dir = '.'.join(os.path.splitext(params['fasta'])[:-1])
  if not os.path.isdir(base_dir):
    os.makedirs(base_dir)

  fasta = "input.fasta"
  shutil.copy(params['fasta'], os.path.join(base_dir, fasta))
  params['fasta'] = fasta

  os.chdir(base_dir)

  # initialize the proteins data structure
  prot_ids, proteins = create_protein_data_structure(fasta)
  features = [signalp4, lipop1, hmmsearch3]
  if 'tmhmm' in params['helix_programs']:
    features.append(tmhmm)
  if 'memsat3' in params['helix_programs']:
    features.append(memsat3)
  for extract_protein_feature in features:
    extract_protein_feature(params, proteins)

  for prot_id in prot_ids:
    details, category = \
        predict_surface_exposure(params, proteins[prot_id])
    if details.endswith(';'):
      details = details[:-1]
    if details is '':
      details = "."
    proteins[prot_id]['details'] = details
    proteins[prot_id]['category'] = category
        
  return prot_ids, proteins


description = """
Inmembrane is a python script that sequentially carries out
bioinformatic analysis of a fasta file, collates the results
and generates a combined analysis of all the analyses.

(c) 2011 Bosco Ho and Andrew Perry
"""

if __name__ == "__main__":
  parser = OptionParser()
  (options, args) = parser.parse_args()

  params = get_params()
  
  if ('fasta' not in params or not params['fasta']) and not args:
      print description
      parser.print_help()
      sys.exit(1)

  if 'fasta' not in params or not params['fasta']:
      params['fasta'] = args[0]

  prot_ids, proteins = identify_pse_proteins(params)

  for prot_id in prot_ids:
    protein = proteins[prot_id]
    word = prot_id.split()[0]
    if "!" in word:
      prot_id = word.split("|")[1]
    else:
      prot_id = word
    print "%-15s %-13s %-50s %s" % \
        (word, 
         protein['category'], 
         protein['details'],
         protein['name'][:60])

