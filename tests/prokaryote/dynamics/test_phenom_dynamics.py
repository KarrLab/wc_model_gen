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

        cls.kb = wc_kb.io.Reader().run('tests/fixtures/min_kb.xlsx',
                                       'tests/fixtures/min_kb_seq.fna',
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

        """ Analysis
        run_results = RunResults(run_results_dir)
        rna_ids=[]
        df = run_results.get('populations')

        for rna in model.species_types.get(type = wc_lang.SpeciesTypeType.rna):
            rna_ids.append(rna.species[0].id())

        txt =  'ATP has depleted to: {} \n'.format(round(df.loc[100.0,'atp[c]']/df.loc[0.0,'atp[c]'],2))
        txt += 'CTP has depleted to: {} \n'.format(round(df.loc[100.0,'ctp[c]']/df.loc[0.0,'ctp[c]'],2))
        txt += 'GTP has depleted to: {} \n'.format(round(df.loc[100.0,'gtp[c]']/df.loc[0.0,'gtp[c]'],2))
        txt += 'UTP has depleted to: {} \n \n'.format(round(df.loc[100.0,'utp[c]']/df.loc[0.0,'utp[c]'],2))
        txt += 'Init copy number mean={}; std={} \n'.format(round(np.mean(df.loc[0.0,rna_ids].values),2),
                                                            round(np.std(df.loc[0.0,rna_ids].values),2))
        txt += 'Final copy number mean={}; std={}'.format(round(np.mean(df.loc[100.0,rna_ids].values),2),
                                                          round(np.std(df.loc[100.0,rna_ids].values),2))
        print(txt)
        """
