import os
import unittest
import sys

# hack to allow tests to find inmembrane in directory above
module_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(module_dir, '..'))

import inmembrane 
import plugins

class TestTmbetadisc_rbf(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(module_dir, 'tatfind')

  def test_tmbetadisc_rbf(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)

    inmembrane.silence_log(True)
    inmembrane.clean_directory('.', ['input.fasta'])
    
    self.params = inmembrane.get_params()
    self.params['fasta'] = "input.fasta"
    self.seqids, self.proteins = \
        inmembrane.create_proteins_dict(self.params['fasta'])
    plugins.tmbetadisc_rbf_web.annotate(self.params, self.proteins)

    self.expected_output = {
      'sp|P31550': {
        'sequence_length': 327, 
        'name': 'THIB_ECOLI Thiamine-binding periplasmic protein OS=Escherichia coli (strain K12) GN=thiB PE=1 SV=2', 
        'seq': 'MLKKCLPLLLLCTAPVFAKPVLTVYTYDSFAADWGPGPVVKKAFEADCNCELKLVALEDGVSLLNRLRMEGKNSKADVVLGLDNNLLDAASKTGLFAKSGVAADAVNVPGGWNNDTFVPFDYGYFAFVYDKNKLKNPPQSLKELVESDQNWRVIYQDPRTSTPGLGLLLWMQKVYGDDAPQAWQKLAKKTVTVTKGWSEAYGLFLKGESDLVLSYTTSPAYHILEEKKDNYAAANFSEGHYLQVEVAARTAASKQPELAQKFLQFMVSPAFQNAIPTGNWMYPVANVTLPAGFEKLTKPATTLEFTPAEVAAQRQAWISEWQRAVSR', 
        'is_tmbetadisc_rbf': False, 
      },
      'sp|P0AB06': {
        'sequence_length': 182, 
        'name': 'YCBK_ECOLI Uncharacterized protein ycbK OS=Escherichia coli (strain K12) GN=ycbK PE=4 SV=1', 
        'seq': 'MDKFDANRRKLLALGGVALGAAILPTPAFATLSTPRPRILTLNNLHTGESIKAEFFDGRGYIQEELAKLNHFFRDYRANKIKSIDPGLFDQLYRLQGLLGTRKPVQLISGYRSIDTNNELRARSRGVAKKSYHTKGQAMDFHIEGIALSNIRKAALSMRAGGVGYYPRSNFVHIDTGPARHW', 
        'is_tmbetadisc_rbf': True, 
      },
      'sp|P31549': {
        'sequence_length': 536, 
        'name': 'THIP_ECOLI Thiamine transport system permease protein thiP OS=Escherichia coli (strain K12) GN=thiP PE=3 SV=2', 
        'seq': 'MATRRQPLIPGWLIPGVSATTLVVAVALAAFLALWWNAPQDDWVAVWQDSYLWHVVRFSFWQAFLSALLSVIPAIFLARALYRRRFPGRLALLRLCAMTLILPVLVAVFGILSVYGRQGWLATLCQSLGLEWTFSPYGLQGILLAHVFFNLPMASRLLLQALENIPGEQRQLAAQLGMRSWHFFRFVEWPWLRRQIPPVAALIFMLCFASFATVLSLGGGPQATTIELAIYQALSYDYDPARAAMLALLQMVCCLGLVLLSQRLSKAIAPGTTLLQGWRDPDDRLHSRICDTVLIVLALLLLLPPLLAVIVDGVNRQLPEVLAQPVLWQALWTSLRIALAAGVLCVVLTMMLLWSSRELRARQKMLAGQVLEMSGMLILAMPGIVLATGFFLLLNNTIGLPQSADGIVIFTNALMAIPYALKVLENPMRDITARYSMLCQSLGIEGWSRLKVVELRALKRPLAQALAFACVLSIGDFGVVALFGNDDFRTLPFYLYQQIGSYRSQDGAVTALILLLLCFLLFTVIEKLPGRNVKTD', 
        'is_tmbetadisc_rbf': False, 
      },
    }

    import helpers
    helpers.print_proteins(self.proteins)
    
    self.assertEqual(self.proteins, self.expected_output)

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
