""" Test of running models with phenomenological rate laws

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-10-01
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults
from wc_model_gen import prokaryote
import numpy as np
import unittest
import wc_lang
import wc_kb_gen
import wc_kb
import unittest
import tempfile
import shutil
import os

class PhenomenologicalTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_growth_transcription(self):

        # Construct model
        self.kb = wc_kb.io.Reader().run('tests/fixtures/min_model_kb.xlsx',
                                       'tests/fixtures/min_model_kb_seq.fna',
                                        strict=False)

        self.model = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = self.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.TranscriptionSubmodelGenerator,
                                              prokaryote.RnaDegradationSubmodelGenerator,
                                              prokaryote.MetabolismSubmodelGenerator],
                        options = {'component': {
                                    'TranscriptionSubmodelGenerator': {
                                        'rate_dynamics': 'phenomenological'},
                                    'RnaDegradationSubmodelGenerator': {
                                        'rate_dynamics': 'phenomenological'}}}).run()

        # Simulate model
        checkpoint_period = 10
        end_time = 100
        simulation = Simulation(self.model)
        results = simulation.run(end_time, self.dir, checkpoint_period)
        self.assertIsInstance(results, tuple)

        # Check results
        num_events  = results[0]
        run_results_dir = results[1]
        run_results = RunResults(run_results_dir)
        df = run_results.get('populations')

        cytosol = self.model.compartments.get_one(id='c')
        rna_ids=[]
        for rna in self.model.species_types.get(type = wc_lang.SpeciesTypeType.rna):
            rna_ids.append(rna.species.get_one(compartment=cytosol).id)

        avg_init_rna_cn  = np.mean(df.loc[0.0,rna_ids].values)
        avg_final_rna_cn = np.mean(df.loc[100.0,rna_ids].values)

        print(simulation.provide_event_counts())
        print('\n INIT AVG RNA COPY NUMBERS:', avg_init_rna_cn)
        print('\n FINAL AVG RNA COPY NUMBERS:', avg_final_rna_cn)

        self.assertTrue(abs(2*avg_init_rna_cn-avg_final_rna_cn) < avg_init_rna_cn*0.15)

class MechanisticDynamicsTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_growth_transcription(self):

        # Construct model
        self.kb = wc_kb.io.Reader().run('tests/fixtures/min_model_kb.xlsx',
                                       'tests/fixtures/min_model_kb_seq.fna',
                                        strict=False)

        self.model = prokaryote.ProkaryoteModelGenerator(
                    knowledge_base = self.kb,
                    component_generators=[prokaryote.InitalizeModel,
                                          prokaryote.TranscriptionSubmodelGenerator,
                                          prokaryote.RnaDegradationSubmodelGenerator,
                                          prokaryote.MetabolismSubmodelGenerator],
                    options = {'component': {
                                'TranscriptionSubmodelGenerator': {
                                    'rate_dynamics': 'mechanistic'},
                                'RnaDegradationSubmodelGenerator': {
                                    'rate_dynamics': 'mechanistic'}}}).run()

        # Simulate model
        checkpoint_period = 10
        end_time = 100
        simulation = Simulation(self.model)
        results = simulation.run(end_time, self.dir, checkpoint_period)
        self.assertIsInstance(results, tuple)

        # Check results
        num_events  = results[0]
        run_results_dir = results[1]
        run_results = RunResults(run_results_dir)
        df = run_results.get('populations')

        cytosol = self.model.compartments.get_one(id='c')
        rna_ids=[]
        for rna in self.model.species_types.get(type = wc_lang.SpeciesTypeType.rna):
            rna_ids.append(rna.species.get_one(compartment=cytosol).id)

        avg_init_rna_cn  = np.mean(df.loc[0.0,rna_ids].values)
        avg_final_rna_cn = np.mean(df.loc[100.0,rna_ids].values)

        print(simulation.provide_event_counts())
        print('\n INIT AVG RNA COPY NUMBERS:', avg_init_rna_cn)
        print('\n FINAL AVG RNA COPY NUMBERS:', avg_final_rna_cn)

        self.assertTrue(abs(2*avg_init_rna_cn-avg_final_rna_cn) < avg_init_rna_cn*0.15)
