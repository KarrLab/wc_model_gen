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

        cls.kb = wc_kb.io.Reader().run('/tests/fixtures/min_model_kb.xlsx',
                                       '/tests/fixtures/min_model_kb_seq.fna',
                                        strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        options = {'component': {
                                    'TranscriptionSubmodelGenerator': {'rate_dynamics': 'phenomenological'},
                                    'RnaDegradationSubmodelGenerator': {'rate_dynamics': 'phenomenological'},
                                    'TranslationSubmodelGenerator': {'rate_dynamics': 'phenomenological'},
                                    'ProteinDegradationSubmodelGenerator': {'rate_dynamics': 'phenomenological'}}}).run()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_exponential_growth(self):
        checkpoint_period = 5
        end_time = 105
        results_dir = self.dir
        simulation = Simulation(self.model)

        results = simulation.run(end_time, results_dir, checkpoint_period)
        num_events  = results[0]
        run_results_dir = results[1]

        self.assertIsInstance(results, tuple)
