import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import tatfind_web

class TestTatfind(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(
       os.path.abspath(
       os.path.dirname(inmembrane.tests.__file__)), 'tatfind')

  def test_tatfind(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)

    helpers.silence_log(True)
    helpers.clean_directory('.', ['input.fasta'])
    
    self.params = inmembrane.get_params()
    self.params['fasta'] = "input.fasta"
    self.seqids, self.proteins = \
        helpers.create_proteins_dict(self.params['fasta'])
    tatfind_web.annotate(self.params, self.proteins)

    self.expected_output = {
      'sp|P31550': {
        'sequence_length': 327, 
        'name': 'THIB_ECOLI Thiamine-binding periplasmic protein OS=Escherichia coli (strain K12) GN=thiB PE=1 SV=2', 
        'seq': 'MLKKCLPLLLLCTAPVFAKPVLTVYTYDSFAADWGPGPVVKKAFEADCNCELKLVALEDGVSLLNRLRMEGKNSKADVVLGLDNNLLDAASKTGLFAKSGVAADAVNVPGGWNNDTFVPFDYGYFAFVYDKNKLKNPPQSLKELVESDQNWRVIYQDPRTSTPGLGLLLWMQKVYGDDAPQAWQKLAKKTVTVTKGWSEAYGLFLKGESDLVLSYTTSPAYHILEEKKDNYAAANFSEGHYLQVEVAARTAASKQPELAQKFLQFMVSPAFQNAIPTGNWMYPVANVTLPAGFEKLTKPATTLEFTPAEVAAQRQAWISEWQRAVSR', 
        'is_tatfind': False, 
      },
      'sp|P0AB06': {
        'sequence_length': 182, 
        'name': 'YCBK_ECOLI Uncharacterized protein ycbK OS=Escherichia coli (strain K12) GN=ycbK PE=4 SV=1', 
        'seq': 'MDKFDANRRKLLALGGVALGAAILPTPAFATLSTPRPRILTLNNLHTGESIKAEFFDGRGYIQEELAKLNHFFRDYRANKIKSIDPGLFDQLYRLQGLLGTRKPVQLISGYRSIDTNNELRARSRGVAKKSYHTKGQAMDFHIEGIALSNIRKAALSMRAGGVGYYPRSNFVHIDTGPARHW', 
        'is_tatfind': True, 
      },
      'sp|P31549': {
        'sequence_length': 536, 
        'name': 'THIP_ECOLI Thiamine transport system permease protein thiP OS=Escherichia coli (strain K12) GN=thiP PE=3 SV=2', 
        'seq': 'MATRRQPLIPGWLIPGVSATTLVVAVALAAFLALWWNAPQDDWVAVWQDSYLWHVVRFSFWQAFLSALLSVIPAIFLARALYRRRFPGRLALLRLCAMTLILPVLVAVFGILSVYGRQGWLATLCQSLGLEWTFSPYGLQGILLAHVFFNLPMASRLLLQALENIPGEQRQLAAQLGMRSWHFFRFVEWPWLRRQIPPVAALIFMLCFASFATVLSLGGGPQATTIELAIYQALSYDYDPARAAMLALLQMVCCLGLVLLSQRLSKAIAPGTTLLQGWRDPDDRLHSRICDTVLIVLALLLLLPPLLAVIVDGVNRQLPEVLAQPVLWQALWTSLRIALAAGVLCVVLTMMLLWSSRELRARQKMLAGQVLEMSGMLILAMPGIVLATGFFLLLNNTIGLPQSADGIVIFTNALMAIPYALKVLENPMRDITARYSMLCQSLGIEGWSRLKVVELRALKRPLAQALAFACVLSIGDFGVVALFGNDDFRTLPFYLYQQIGSYRSQDGAVTALILLLLCFLLFTVIEKLPGRNVKTD', 
        'is_tatfind': True, 
      },
    }

    #import helpers
    #helpers.print_proteins(self.proteins)
    
    self.assertEqual(self.proteins, self.expected_output)

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
