import glob
import math
import os
import sys

parms = {
  'organism': 'gram+',
  'signalp4_bin': 'signalp',
  'lipop1_bin': 'LipoP',
  'tmhmm_bin': 'tmhmm',
  'hmmsearch3_bin': 'hmmsearch',
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
          profile_name = os.path.basename(parms['hmm_profile'])
          profile_name = profile_name.replace('.hmmsearch', '')
          proteins[name]['hmmsearch'].append(profile_name)


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
    return "PSE, %s profile match" % protein['hmmsearch']

  if protein['is_lipop']: 
    chop_nterminal_peptide(protein, protein['lipop_cleave_position'])
    if protein['n_tmhmm_helix'] == 0:
      if protein['sequence_length'] < parms['terminal_exposed_loop_min']:
        return "MEMBRANE"
      else:
        return "PSE"
  elif protein['is_signalp']:
    chop_nterminal_peptide(protein, protein['signalp_cleave_position'])
    if protein['n_tmhmm_helix'] == 0:
      return "SECRETED"
  if protein['n_tmhmm_helix'] > 0:
    if has_surface_exposed_loop(parms, protein):
      return "PSE"
    else:
      return "MEMBRANE"
  return "CYTOPLASMIC"


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
  for extract_protein_feature in \
      [signalp4, lipop1, tmhmm, hmmsearch3]:
    extract_protein_feature(parms, proteins)
  for prot_id in prot_ids:
    proteins[prot_id]['category'] = \
        predict_surface_exposure(parms, proteins[prot_id])
  return prot_ids, proteins



if __name__ == "__main__":
  fasta = sys.argv[1]
  prot_ids, proteins = identify_pse_proteins(
      parms, fasta)
  for prot_id in prot_ids:
    protein = proteins[prot_id]
    print "%-7s %-13s %s" % \
        (prot_id.split("|")[1], 
         protein['category'].split()[0], 
         protein['name'][:60])


  





