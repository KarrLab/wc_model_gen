""" Tests of metabolism submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb_gen import random
from wc_model_gen.prokaryote import metabolism
import unittest
import wc_lang


class MetabolismSubmodelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = random.RandomKbGenerator(options={
            'component': {
                'GenomeGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 200,
                    'mean_gene_len': 10,
                },
                'PropertiesGenerator': {
                    'mean_cell_density': 1e6,
                    'mean_fraction_dry_weight': 0.3
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

        extra = model.compartments.get_one(id='e')
        self.assertEqual(extra.name, 'extracellular space')

        # check parameters generated
        self.assertEqual(model.parameters.get_one(
            id='fractionDryWeight').value, 0.3)

        # check species types and species generated
        atp = model.species_types.get_one(id='atp')
        atp_cytosol = atp.species.get_one(compartment=cytosol)
        self.assertEqual(atp_cytosol.concentration.units,
                         wc_lang.ConcentrationUnit.M)
