""" Tests of metabolism submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import rand_wc_model_gen
from rand_wc_model_gen import kb_gen
from wc_model_gen.rand_gen import metabolism
import unittest
import wc_kb
import wc_lang


class MetabolismSubmodelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = kb_gen.KbGenerator(options={
            'component': {
                'GenomeGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 100,
                    'mean_gene_len': 100,
                },
                'PropertiesGenerator': {
                    'mean_cell_density': 1e6,
                    'mean_fraction_dry_weight': 0.3,
                },
            },
        }).run()
        cell = kb.cell

        model = wc_lang.Model()
        gen = metabolism.MetabolismSubmodelGenerator(kb, model, options={})
        gen.run()

        submodel = model.submodels.get_one(id='metabolism')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')

        cytosol = model.compartments.get_one(id='e')
        self.assertEqual(cytosol.name, 'extracellular space')

        # check parameters generated
        self.assertEqual(model.parameters.get_one(id='fractionDryWeight').value, 0.3)
