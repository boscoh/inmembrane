import glob
import math
import os
import sys

parms = {
  'organism': 'gram+',
  'signalp4_bin': '/home/boscoh/packages/signalp-4.0/signalp',
  'lipop1_bin': 'LipoP',
  'tmhmm_bin': 'tmhmm',
  'hmmsearch3_bin': '/bio/sw/hmmer-3.0-linux-intel-x86_64/bin/hmmsearch',
  'hmm_profiles_dir': 'hmm_profiles',
  'hmm_evalue_cutoff': 0.01,
  'terminal_exposed_loop_min': 50,
  'internal_exposed_loop_min': 100,
}


def basename(parms):
  return '.'.join(os.path.splitext(parms['fasta'])[:-1])


def run(cmd, out_file="log.txt"):
  full_cmd = cmd + " > " + out_file
  print "#", full_cmd
  if os.path.isfile(out_file):
    print "# -> skipped: %s already exists" % out_file
  else:
    os.system(full_cmd)


def read_fasta_keys(fname):
  return [l[1:] for l in open(fname) if l.startswith(">")]


def hmmsearch3(parms, proteins):
  base = basename(parms)
  profiles_dir = parms['hmm_profiles_dir']

  for hmm_profile in glob.glob(profiles_dir + '/*.hmm'):
    parms['hmm_profile'] = hmm_profile
    hmm_profile = os.path.basename(parms['hmm_profile'])
    hmm_name = hmm_profile.replace('.hmm', '')
    hmmsearch3_out = '%s.hmm.%s.out' % (base, hmm_name)
    run('%(hmmsearch3_bin)s -Z 2000 -E 10 %(hmm_profile)s %(fasta)s' % parms, hmmsearch3_out)
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

def has_surface_exposed_loop(parms, protein):
  terminal_exposed_loop_min = parms['terminal_exposed_loop_min']
  internal_exposed_loop_min = parms['internal_exposed_loop_min']

  tmhmm_outer_loops = protein['tmhmm_outer_loops']
  sequence_length = protein['sequence_length']
  has_no_transmembrane_helices = protein['n_tmhmm_helix'] == 0

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
    if has_surface_exposed_loop(parms, protein):
      return s, "PSE"
    else:
      return s, "MEMBRANE"

  return s, "CYTOPLASM"

def init(parms):
  """
  Initialize the proteins data structure (dictionary keyed by sequence id).
  Takes a parameter dictionary as input.
  Returns a tuple (sequence id, proteins).
  """
  # initialize the proteins data structure
  headers = read_fasta_keys(parms['fasta'])
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
    
  return prot_ids, proteins

def identify_pse_proteins(parms):
  prot_ids, proteins = init(parms)
  
  for extract_protein_feature in \
      [signalp4, lipop1, tmhmm, hmmsearch3, bomp_web]:
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
  parms['fasta'] = fasta
  prot_ids, proteins = identify_pse_proteins(parms)
  for prot_id in prot_ids:
    protein = proteins[prot_id]
    word = prot_id.split()[0]
    if "!" in word:
      prot_id = word.split("|")[1]
    else:
      prot_id = word
    print "%-15s %-13s %-20s %s" % \
        (word, 
         protein['category'], 
         protein['details'],
         protein['name'][:60])

