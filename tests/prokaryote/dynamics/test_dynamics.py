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

class PhenomDynamicsTestCase(unittest.TestCase):

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
                                              prokaryote.MetabolismSubmodelGenerator],
                        options = {'component': {
                                    'TranscriptionSubmodelGenerator': {
                                        'rate_dynamics': 'phenomenological'},
                                    'RnaDegradationSubmodelGenerator': {
                                        'rate_dynamics': 'phenomenological'}}}).run()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_exponential_growth(self):
        checkpoint_period = 10
        end_time = 100

        simulation = Simulation(self.model)
        results = simulation.run(end_time, self.dir, checkpoint_period)
        self.assertIsInstance(results, tuple)

        num_events  = results[0]
        run_results_dir = results[1]

        run_results = RunResults(run_results_dir)
        rna_ids=[]
        df = run_results.get('populations')
        for rna in self.model.species_types.get(type = wc_lang.SpeciesTypeType.rna):
            rna_ids.append(rna.species[0].id)

        avg_init_rna_cn  = np.mean(df.loc[0.0,rna_ids].values)
        avg_final_rna_cn = np.mean(df.loc[100.0,rna_ids].values)

        #print(simulation.provide_event_counts())
        #print('\n INIT:', avg_init_rna_cn)
        #print('\n FINAL:', avg_final_rna_cn)

        # Check if RNA content has doubled after CC - 15% tolerance
        self.assertTrue(abs(2*avg_init_rna_cn-avg_final_rna_cn) < avg_init_rna_cn*0.2)

class MechanisticDynamicsTestCase(unittest.TestCase):

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
                                          prokaryote.MetabolismSubmodelGenerator],
                    options = {'component': {
                                'TranscriptionSubmodelGenerator': {
                                    'rate_dynamics': 'mechanistic'},
                                'RnaDegradationSubmodelGenerator': {
                                    'rate_dynamics': 'mechanistic'}}}).run()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_exponential_growth(self):
        checkpoint_period = 10
        end_time = 100

        simulation = Simulation(self.model)
        results = simulation.run(end_time, self.dir, checkpoint_period)
        self.assertIsInstance(results, tuple)

        num_events  = results[0]
        run_results_dir = results[1]

        run_results = RunResults(run_results_dir)
        rna_ids=[]
        df = run_results.get('populations')
        for rna in self.model.species_types.get(type = wc_lang.SpeciesTypeType.rna):
            rna_ids.append(rna.species[0].id)

        avg_init_rna_cn  = np.mean(df.loc[0.0,rna_ids].values)
        avg_final_rna_cn = np.mean(df.loc[100.0,rna_ids].values)

        #print(simulation.provide_event_counts())
        #print('\n INIT:', avg_init_rna_cn)
        #print('\n FINAL:', avg_final_rna_cn)

        # Check if RNA content has doubled after CC - 15% tolerance
        self.assertTrue(abs(2*avg_init_rna_cn-avg_final_rna_cn) < avg_init_rna_cn*0.2)
