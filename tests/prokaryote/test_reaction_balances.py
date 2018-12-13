""" Test element and charge balance of reactions

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-12-05
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_model_gen import prokaryote
import unittest
import wc_kb
import wc_lang


class ReactionBalanceTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/min_model_kb.xlsx',
                                       'tests/fixtures/min_model_kb_seq.fna',
                                       strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=cls.kb,
            component_generators=[prokaryote.InitalizeModel,
                                  prokaryote.TranscriptionSubmodelGenerator,
                                  prokaryote.RnaDegradationSubmodelGenerator,
                                  prokaryote.MetabolismSubmodelGenerator,
                                  ]).run()

    def test(self):
        self.assertEqual(wc_lang.Validator().run(self.model), None)
