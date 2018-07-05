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

class DegradationSubmodelTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._results_dir = tempfile.mkdtemp()

        cls.kb = wc_kb.io.Reader().run('/home/balazs/Desktop/mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/core.xlsx',
        '/home/balazs/Desktop/mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/seq.fna', strict = False)

        cls.model = wc_lang.io.Reader().run('/home/balazs/Desktop/wc_test/tests/fixtures/min_model.xlsx')

        # There is sth wrong with the output format, on todo list
        #cls.model = wc_model_gen.ModelGenerator(knowledge_base=cls.kb,
        #                                        component_generators=[prokaryote.CompartmentsGenerator,
        #                                                              prokaryote.ParametersGenerator,
        #                                                              prokaryote.MetaboliteSpeciesGenerator,
        #                                                              prokaryote.DegradationSubmodelGenerator]).run()

        cls.degradation_static_test_case  = wc_test.StaticTestCase(model=cls.model)
        cls.degradation_dynamic_test_case = wc_test.DynamicTestCase(model=cls.model)
        cls.run_results = cls.degradation_dynamic_test_case.simulate(end_time=300)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._results_dir)

    def test_expected_classes(self):
        self.assertIsInstance(self.degradation_static_test_case, wc_test.StaticTestCase)
        self.assertIsInstance(self.degradation_dynamic_test_case, wc_test.DynamicTestCase)
        self.assertIsInstance(self.run_results, list)
        self.assertIsInstance(self.run_results[0], wc_sim.multialgorithm.run_results.RunResults)

    def test_concentrations(self):
        pass

    def test_rate_law_scan(self):
        pass
