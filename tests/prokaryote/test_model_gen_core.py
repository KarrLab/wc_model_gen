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
import wc_lang
import wc_utils.util.string


class ModelGeneratorTestCase(unittest.TestCase):
    def test(self):

        rand_kb = random.RandomKbGenerator(options={
                     'component': {
                         'GenomeGenerator': {
                             'mean_num_genes': 30,
                             'mean_gene_len': 50,
                             'num_ncRNA': 0,
                             'translation_table': 4,
                             'mean_copy_number': 200,
                             'mean_half_life': 100},
                         'PropertiesGenerator': {
                             'mean_cell_cycle_length': 100,},
                     },
                 }).run()

        model = prokaryote.ProkaryoteModelGenerator(rand_kb).run()
        self.assertIsInstance(model.submodels.get_one(id='metabolism'), wc_lang.Submodel)

        errors = obj_model.Validator().run(model, get_related=True)
        self.assertEqual(errors, None, msg=wc_utils.util.string.indent_forest(errors))
