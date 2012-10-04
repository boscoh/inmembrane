# -*- coding: utf-8 -*-
from inmembrane.helpers import run, parse_fasta_header, dict_get

citation = {'ref': u"Agnieszka S. Juncker, Hanni Willenbrock, "
                   u"Gunnar Von Heijne, Søren Brunak, Henrik Nielsen, "
                   u"And Anders Krogh. (2003) Prediction of lipoprotein "
                   u"signal peptides in Gram-negative bacteria. Protein "
                   u"Science 12:1652–1662. \n"
                   u"<http://dx.doi.org/10.1110/ps.0303703>",
            'name': 'LipoP 1.0'
           } 


def annotate(params, proteins):
  """
  Uses LipoP to identify lipo-attachment signals in the protein. 
  The 'proteins' dictionary is annotated by adding two fields:
    - 'is_lipop' is a boolean indicating whether a signal is found or not
    - 'lipop_cleave_position' gives the position where the signal
        is cleaved and the protein is attached to a lipid

  Returns a reference to the proteins data structure.
  """

  lipop1_out = 'lipop.out'
  run('%(lipop1_bin)s %(fasta)s' % params, lipop1_out)

  # initialize fields in each protein
  for seqid in proteins:
    proteins[seqid]['is_lipop'] = False
    proteins[seqid]['lipop_cleave_position'] = None

  for l in open(lipop1_out):
    words = l.split()

    if 'SpII score' in l:
      seqid = parse_fasta_header(words[1])[0]
      if 'cleavage' in l:
        pair = words[5].split("=")[1]
        i = int(pair.split('-')[0])
      else:
        i = None
      proteins[seqid]['is_lipop'] = 'Sp' in words[2]
      proteins[seqid]['lipop_cleave_position'] = i

    # check for an E.coli style inner membrane retention signal
    # Asp+2 to cleavage site. There are other apparent retention 
    # signals in E. coli and other gram- bacteria in addition to
    # the Asp+2 which we don't detect here (yet).
    # (Yamaguchi et al, 1988; Tokuda and Matsuyama, 2005 [review])
    if dict_get(proteins[seqid], 'lipop_cleave_position'):
      plus_two = proteins[seqid]['lipop_cleave_position']+1
      if proteins[seqid]['seq'][plus_two] == 'D':
        proteins[seqid]['lipop_im_retention_signal'] = True
      
  return proteins
