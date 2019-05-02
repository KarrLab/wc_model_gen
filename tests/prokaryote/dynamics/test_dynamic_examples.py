""" Testing that toy models behave as expected

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2019-01-25
:Copyright: 2019, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
from wc_sim.multialgorithm.run_results import RunResults
from wc_sim.multialgorithm.simulation import Simulation
from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import numpy as np
import shutil
import tempfile
import unittest
import wc_kb_gen
import wc_lang
import wc_kb
import os


class DynamicsTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp_dirname = tempfile.mkdtemp()

        cls.model = wc_lang.core.Model(id='flat_rate_transfer_model', version='1.0.0')
        model = cls.model

        # Create submodels
        submodel = model.submodels.create(id='transfer')

        # Construct compartments
        comp_A = model.compartments.create(id='A', name='A', mean_init_volume=1E-15)
        comp_A.init_density = model.parameters.create(id='density_A', value=1100., units=unit_registry.parse_units('g l^-1'))
        comp_B = model.compartments.create(id='B', name='B', mean_init_volume=1E-15)
        comp_B.init_density = model.parameters.create(id='density_B', value=1100., units=unit_registry.parse_units('g l^-1'))

        # Construct species / species types
        species_type1 = model.species_types.create(id='st1', charge = 0)
        specie1A = model.species.create(species_type=species_type1, compartment=comp_A)
        specie1B = model.species.create(species_type=species_type1, compartment=comp_B)
        specie1A.id = specie1A.gen_id()
        specie1B.id = specie1B.gen_id()

        species_type2 = model.species_types.create(id='st2', charge = 0)
        specie2A = model.species.create(species_type=species_type2, compartment=comp_A)
        specie2B = model.species.create(species_type=species_type2, compartment=comp_B)
        specie2A.id = specie2A.gen_id()
        specie2B.id = specie2B.gen_id()

        species_type3 = model.species_types.create(id='st3', charge = 0)
        specie3A = model.species.create(species_type=species_type3, compartment=comp_A)
        specie3B = model.species.create(species_type=species_type3, compartment=comp_B)
        specie3A.id = specie3A.gen_id()
        specie3B.id = specie3B.gen_id()

        # Define init concentration distributions
        model.distribution_init_concentrations.create(
            id = 'dist-init-conc-st1[A]',
            species = specie1A,
            mean = 300,
            std = 0,
            units = unit_registry.parse_units('molecules'))
        model.distribution_init_concentrations.create(
            id = 'dist-init-conc-st1[B]',
            species = specie1B,
            mean = 0,
            std = 0,
            units = unit_registry.parse_units('molecules'))
        model.distribution_init_concentrations.create(
            id = 'dist-init-conc-st2[A]',
            species = specie2A,
            mean = 0,
            std = 0,
            units = unit_registry.parse_units('molecules'))
        model.distribution_init_concentrations.create(
            id = 'dist-init-conc-st2[B]',
            species = specie2B,
            mean = 200,
            std = 0,
            units = unit_registry.parse_units('molecules'))
        model.distribution_init_concentrations.create(
            id = 'dist-init-conc-st3[A]',
            species = specie3A,
            mean = 100,
            std = 0,
            units = unit_registry.parse_units('molecules'))
        model.distribution_init_concentrations.create(
            id = 'dist-init-conc-st3[B]',
            species = specie3B,
            mean = 100,
            std = 0,
            units = unit_registry.parse_units('molecules'))

        # Define reactions
        reaction1 = model.reactions.create(submodel=submodel, id='transfer_1A1B', name = 'transfer_1A1B', participants =[])
        reaction1.participants.add(specie1A.species_coefficients.get_or_create(coefficient=-1))
        reaction1.participants.add(specie1B.species_coefficients.get_or_create(coefficient=1))

        reaction2 = model.reactions.create(submodel=submodel, id='transfer_2B2A', name = 'transfer_2B2A', participants =[])
        reaction2.participants.add(specie2A.species_coefficients.get_or_create(coefficient=1))
        reaction2.participants.add(specie2B.species_coefficients.get_or_create(coefficient=-1))

        reaction3 = model.reactions.create(submodel=submodel, id='transfer_3A3B', name = 'transfer_3A3B', participants =[], reversible = True)
        reaction3.participants.add(specie3A.species_coefficients.get_or_create(coefficient=1))
        reaction3.participants.add(specie3B.species_coefficients.get_or_create(coefficient=-1))

        # Define rate laws
        k_cat_r1 = model.parameters.create(id='k_cat_r1', value=1, units=unit_registry.parse_units('1 / second'))
        k_cat_r2 = model.parameters.create(id='k_cat_r2', value=2, units=unit_registry.parse_units('1 / second'))
        k_cat_r3 = model.parameters.create(id='k_cat_r3', value=4, units=unit_registry.parse_units('1 / second'))

        rate_law1 = model.rate_laws.create(id='transfer_1A1B-forward', reaction=reaction1, direction=wc_lang.RateLawDirection.forward)
        rate_law2 = model.rate_laws.create(id='transfer_2B2A-forward', reaction=reaction2, direction=wc_lang.RateLawDirection.forward)
        rate_law3f = model.rate_laws.create(id='transfer_3A3B-forward', reaction=reaction3, direction=wc_lang.RateLawDirection.forward)
        rate_law3b = model.rate_laws.create(id='transfer_3A3B-backward', reaction=reaction3, direction=wc_lang.RateLawDirection.backward)

        objects = {wc_lang.Parameter: {
                         k_cat_r1.id: k_cat_r1,
                         k_cat_r2.id: k_cat_r2,
                         k_cat_r3.id: k_cat_r3},
                   wc_lang.Species: {
                        specie1A.id: specie1A,
                        specie1B.id: specie1B,
                        specie2A.id: specie2A,
                        specie2B.id: specie2B,
                        specie3A.id: specie3A,
                        specie3B.id: specie3B}}

        rate_law1.expression, error = wc_lang.RateLawExpression.deserialize('k_cat_r1', objects)
        rate_law2.expression, error = wc_lang.RateLawExpression.deserialize('k_cat_r2', objects)
        rate_law3f.expression, error = wc_lang.RateLawExpression.deserialize('k_cat_r3', objects)
        rate_law3b.expression, error = wc_lang.RateLawExpression.deserialize('k_cat_r3', objects)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_dirname)

    def test_flat_rates(self):
        model = self.model

        end_time          = 400
        checkpoint_period = 10
        simulation = Simulation(model)

        # Run the simulation 20 times
        #run_results_dir = []
        #for i in range(0,20):
        results = simulation.run(end_time, self.tmp_dirname, checkpoint_period)

        self.assertIsInstance(results, tuple)
        self.assertTrue(len(results)==2)
