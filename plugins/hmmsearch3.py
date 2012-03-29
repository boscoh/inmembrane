import os
import glob
import helpers


def annotate_hmmsearch3(params, proteins):
  """
  Returns a reference to the proteins data structure.

  Uses HMMER to identify sequence motifs in proteins. This function
  annotates the proteins with:
    - 'hmmsearch': a list of motifs that are found in the protein. The
       motifs correspond to the basename of the .hmm files found in the directory
       indicated by the 'hmm_profiles_dir' field of 'params'.
  """

  helpers.log_stderr("# Searching for HMMER profiles in " + params['hmm_profiles_dir'])

  file_tag = os.path.join(params['hmm_profiles_dir'], '*.hmm')
  for hmm_profile in glob.glob(file_tag):
    params['hmm_profile'] = hmm_profile

    hmm_profile = os.path.basename(params['hmm_profile'])
    hmm_name = hmm_profile.replace('.hmm', '')
    hmmsearch3_out = 'hmm.%s.out' % hmm_name

    cmd = '%(hmmsearch3_bin)s -Z 2000 -E 10 %(hmm_profile)s %(fasta)s' % params
    helpers.run(cmd, hmmsearch3_out)

    # init proteins data structure with blank hmmsearch field first
    for seqid in proteins:
      if 'hmmsearch' not in proteins[seqid]:
        proteins[seqid]['hmmsearch'] = []

    # parse the hmmsearch output file
    seqid = None
    for l in open(hmmsearch3_out):
      words = l.split()

      if l.startswith(">>"):
        seqid = helpers.parse_fasta_header(l[3:])[0]
        continue

      if seqid is None:
        continue

      if 'conditional E-value' in l:
        evalue = float(words[-1])
        score = float(words[-5])
        if evalue <= params['hmm_evalue_max'] and \
            score >= params['hmm_score_min']:
          proteins[seqid]['hmmsearch'].append(hmm_name)


  return proteins
