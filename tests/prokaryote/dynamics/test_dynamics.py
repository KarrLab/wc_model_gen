""" Test of running models with phenomenological rate laws

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-10-01
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
from wc_sim.multialgorithm.run_results import RunResults
from wc_sim.multialgorithm.simulation import Simulation
from wc_onto import onto as wc_ontology
import numpy as np
import os
import shutil
import tempfile
import unittest
import wc_kb
import wc_kb_gen
import wc_lang


class DynamicsTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    @unittest.skip('fix charge/element imbalance')
    def test_growth_transcription(self):

        # Construct model
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_kb__DOT__io__DOT__strict', '0')
        with env:
            self.kb = wc_kb.io.Reader().run('tests/fixtures/min_model_kb.xlsx',
                                            'tests/fixtures/min_model_kb_seq.fna',
                                            )[wc_kb.KnowledgeBase][0]

        self.model = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=self.kb,
            component_generators=[prokaryote.InitalizeModel,
                                  prokaryote.TranscriptionSubmodelGenerator,
                                  prokaryote.RnaDegradationSubmodelGenerator,
                                  prokaryote.MetabolismSubmodelGenerator],
            options={'component': {
                'TranscriptionSubmodelGenerator': {
                    'beta': 1.},
                'RnaDegradationSubmodelGenerator': {
                    'beta': 1.}}}).run()

        # Simulate model
        checkpoint_period = 10
        end_time = 100
        simulation = Simulation(self.model)
        results = simulation.run(end_time, self.dir, checkpoint_period)
        self.assertIsInstance(results, tuple)

        # Check results
        num_events = results[0]
        run_results_dir = results[1]
        run_results = RunResults(run_results_dir)
        df = run_results.get('populations')

        cytosol = self.model.compartments.get_one(id='c')
        rna_ids = []
        for rna in self.model.species_types.get(type=wc_ontology['WC:RNA']): # RNA
            rna_ids.append(rna.species.get_one(compartment=cytosol).id)

        avg_init_rna_cn = np.mean(df.loc[0.0, rna_ids].values)
        avg_final_rna_cn = np.mean(df.loc[100.0, rna_ids].values)

        print(simulation.provide_event_counts())
        print('\n INIT AVG RNA COPY NUMBERS:', avg_init_rna_cn)
        print('\n FINAL AVG RNA COPY NUMBERS:', avg_final_rna_cn)

        self.assertGreater(avg_final_rna_cn, 2 * avg_init_rna_cn - avg_init_rna_cn * 0.15)
        self.assertLess(avg_final_rna_cn, 2 * avg_init_rna_cn + avg_init_rna_cn * 0.15)
