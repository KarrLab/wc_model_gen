""" Tests of model generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""
from wc_kb_gen import random
from wc_model_gen import rand_gen
import obj_model
import unittest
import wc_lang
import wc_utils.util.string


class ModelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = random.RandomKbGenerator(options={
            'component': {
                'GenomeGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 100,
                    'mean_gene_len': 10,
                },
            },
        }).run()
        model = rand_gen.RandomModelGenerator(kb).run()

        self.assertIsInstance(model.submodels.get_one(
            id='transcription'), wc_lang.Submodel)

        errors = obj_model.Validator().run(model, get_related=True)
        self.assertEqual(
            errors, None, msg=wc_utils.util.string.indent_forest(errors))
