""" Test of running models with phenomenological rate laws

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-10-01
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults
import wc_model_gen.prokaryote as prokaryote
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
                                              prokaryote.MetabolismSubmodelGenerator]).run()

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
        for rna in model.species_types.get(type = wc_lang.SpeciesTypeType.rna):
            rna_ids.append(rna.species[0].id())

        avg_init_rna_cn  = np.mean(df.loc[0.0,rna_ids].values)
        avg_final_rna_cn = np.mean(df.loc[100.0,rna_ids].values)

        # Check if RNA content has doubled after CC - 15% tolerance
        self.assertTrue(abs(2*avg_init_rna_cn-avg_final_rna_cn) < avg_init_rna_cn*0.15)
