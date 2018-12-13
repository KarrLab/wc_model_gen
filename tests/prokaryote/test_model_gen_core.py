""" Tests of model generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""
from wc_kb_gen import random
from wc_model_gen import prokaryote
import obj_model
import unittest
import wc_kb
import wc_lang
import wc_utils.util.string


class ModelGeneratorTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/test_broken_kb.xlsx',
                                       'tests/fixtures/test_broken_seq.fna',
                                       strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(knowledge_base=cls.kb).run()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_submodels(self):
        self.assertEqual(5, len(self.model.submodels))

    def test_compartments(self):
        cytosol = self.model.compartments.get_one(id='c')
        extracellular_space = self.model.compartments.get_one(id='e')
        self.assertIsInstance(cytosol, wc_lang.Compartment)
        self.assertIsInstance(extracellular_space, wc_lang.Compartment)

    def test_parameters(self):
        cell_cycle_len = self.model.parameters.get_one(id='cell_cycle_len')
        fraction_dry_weight = self.model.parameters.get_one(id='fraction_dry_weight')
        self.assertIsInstance(cell_cycle_len, wc_lang.Parameter)
        self.assertIsInstance(fraction_dry_weight, wc_lang.Parameter)

        errors = obj_model.Validator().run(self.model, get_related=True)
        self.assertEqual(errors, None, msg=wc_utils.util.string.indent_forest(errors))
