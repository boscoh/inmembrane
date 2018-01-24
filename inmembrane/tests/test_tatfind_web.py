import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import tatfind_web

from inmembrane.tests.PluginTestBase import PluginTestBase


class TestTatfind(PluginTestBase):
    _plugin_name = "tatfind_web"

    def test_tatfind_web(self):
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

        # import helpers
        # helpers.print_proteins(self.proteins)
        for seqid in self.proteins:
            self.assertEqual(self.proteins[seqid]["is_tatfind"],
                             self.expected_output[seqid]["is_tatfind"])


if __name__ == '__main__':
    unittest.main()
