""" Tests of the core of eukaryote model generator

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-09
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import core
from wc_model_gen.eukaryote import initialize_model
import unittest
import wc_kb


class TestCase(unittest.TestCase):

    def test_generator(self):

        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()

        cell.properties.create(id='cell_volume', value=1E-15)
        cell.properties.create(id='mean_doubling_time', value=20., units='hours')

        compartments = {'c': 0.7, 'n': 0.3, 'e': None}
        for k, v in compartments.items():
            cell.compartments.create(id=k, volumetric_fraction=v)

        test_component_generators = [
            initialize_model.InitializeModel
            ]

        gen = core.EukaryoteModelGenerator(kb, component_generators=test_component_generators,
            options={
                'id': 'h1_test',
                'name': 'h1 model',
                'version': '2.3.1',
                'component':{
                    'InitializeModel': {
                        'cell_density': 1100.,
                        'gen_dna': False,
                        'gen_pre_rnas': False,
                        'gen_transcripts': False,
                        'gen_protein': False,
                        'gen_metabolites': False,
                        'gen_complexes': False,
                        'gen_distribution_init_concentrations': False,
                        'gen_observables': False,
                        'gen_kb_reactions': False,
                        'gen_kb_rate_laws': False,
                        }            
                    }	       
                }
            )

        model = gen.run()

        self.assertEqual(model.id, 'h1_test')
        self.assertEqual(model.name, 'h1 model')
        self.assertEqual(model.version, '2.3.1')
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').value, 72000)
        self.assertEqual(model.parameters.get_one(id='density_n').value, 1100)
        self.assertEqual(model.parameters.get_one(id='density_e').value, 1000)
        self.assertEqual(model.compartments.get_one(id='n').mean_init_volume, 3E-16)
        self.assertEqual(model.compartments.get_one(id='e').mean_init_volume, 1E-10)
