import os
import glob
import helpers


def hmmsearch3(params, proteins):
  file_tag = os.path.join(params['hmm_profiles_dir'], '*.hmm')
  helpers.log_stderr("# Searching for HMMER profiles in " + params['hmm_profiles_dir'])
  for hmm_profile in glob.glob(file_tag):
    params['hmm_profile'] = hmm_profile
    hmm_profile = os.path.basename(params['hmm_profile'])
    hmm_name = hmm_profile.replace('.hmm', '')
    hmmsearch3_out = 'hmm.%s.out' % hmm_name
    helpers.run(
        '%(hmmsearch3_bin)s -Z 2000 -E 10 %(hmm_profile)s %(fasta)s' % \
        params, hmmsearch3_out)
    name = None
    for l in open(hmmsearch3_out):
      words = l.split()
      if l.startswith(">>"):
        name = helpers.parse_fasta_header(l[3:])[0]
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
          proteins[name]['hmmsearch'].append(hmm_name)
  return proteins
