""" Tests of transcription submodel generation

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-12-05
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_test import StaticTestCase
from wc_model_gen import prokaryote
import unittest
import wc_kb
import wc_lang
import tempfile
import shutil
import os

class TestReactionBalances(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

        cls.kb = wc_kb.io.Reader().run('tests/fixtures/min_model_kb.xlsx',
                                       'tests/fixtures/min_model_kb_seq.fna',
                                        strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(
                    knowledge_base = cls.kb,
                    component_generators=[prokaryote.InitalizeModel,
                                          prokaryote.TranscriptionSubmodelGenerator,
                                          prokaryote.RnaDegradationSubmodelGenerator,
                                          prokaryote.MetabolismSubmodelGenerator]).run()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    @unittest.skip('Test')
    def test_reaction(self):
        mass_balanced = StaticTestCase(model=self.model).reactions_mass_balanced()
        charge_balanced = StaticTestCase(model=self.model).reactions_charge_balanced()

        print(mass_balanced)
        print(all(mass_balanced))

        self.assertTrue(all(mass_balanced))
        self.assertTrue(all(charge_balanced))
