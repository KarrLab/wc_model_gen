import os
import shutil
import unittest
import tempfile
import wc_kb
import wc_sim
import wc_test
import wc_lang
import wc_model_gen
import wc_model_gen.prokaryote as prokaryote

class TranscriptionSubmodelTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('\n transcription setup CLASS run')
        cls._results_dir = tempfile.mkdtemp()

        cls.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx','tests/fixtures/seq.fna', strict = False)
        cls.model = wc_model_gen.ModelGenerator(
                            knowledge_base=cls.kb,
                            component_generators=[
                               prokaryote.CompartmentsGenerator,
                               prokaryote.ParametersGenerator,
                               prokaryote.MetaboliteSpeciesGenerator,
                               prokaryote.TranscriptionSubmodelGenerator]).run()

        cls.model.id = 'prokaryote_submodels_tests'
        cls.model.version = '0.0.1'

    @classmethod
    def tearDownClass(cls):
        print('\n transcription teardown CLASS run')
        shutil.rmtree(cls._results_dir)

    def tearDown(cls):
        print('\n transcription teardown run')

    def test_submodel_added(self):
        submodel_ids = []
        for submodel in self.model.get_submodels():
            submodel_ids.append(submodel.id)

        self.assertTrue('transcription' in submodel_ids)

    # Check if transcription proceeds without triphosphate molecules
    def test_are_triphosphates_bottleneck_to_transcription(self):

        mod_species={'ATP[c]':0,'CTP[c]':0,'GTP[c]':0,'UTP[c]':0}

        rna_species=[]
        for specie_type in self.model.species_types.get(type=4):
            for specie in specie_type.species:
                rna_species.append(specie.id())

        test_case = wc_test.DynamicTestCase(model=self.model).perturb_species(mod_species=mod_species)
        run_results = test_case.simulate(end_time=900)

        for rna_specie in rna_species:
            concentration = run_results.get('populations')[rna_specie].values

            self.assertEqual(concentration[len(concentration)-1],concentration[0])



"""
class TranslationSubmodelTestCase(ProkaryoteSubmodelsTestCases):

    def setUp(self):
        print('\n SETUP RUN')
        self.model = prokaryote.TranslationSubmodelGenerator(self.kb, self.base_model).run()

class DegradationSubmodelTestCase(ProkaryoteSubmodelsTestCases):

    def setUp(self):
        print('\n SETUP RUN')
        self.model = prokaryote.DegradationSubmodelGenerator(self.kb, self.base_model).run()
"""

"""
metabolite = 1
protein = 2
dna = 3
rna = 4
pseudo_species = 5

self.degradation_static_test_case  = wc_test.StaticTestCase(model=self.base_model)
self.degradation_dynamic_test_case = wc_test.DynamicTestCase(model=self.base_model)
self.run_results = self.degradation_dynamic_test_case.simulate(end_time=300)

# Test expected classes
self.assertIsInstance(self.degradation_static_test_case, wc_test.StaticTestCase)
self.assertIsInstance(self.degradation_dynamic_test_case, wc_test.DynamicTestCase)
self.assertIsInstance(self.run_results, list)
self.assertIsInstance(self.run_results[0], wc_sim.multialgorithm.run_results.RunResults)
"""
