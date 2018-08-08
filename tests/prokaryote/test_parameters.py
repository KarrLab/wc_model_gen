import shutil
import os
import unittest
import tempfile
import wc_kb
import wc_lang
import wc_model_gen
import wc_model_gen.prokaryote as prokaryote


class ParametersGeneratorTestCase(unittest.TestCase):

    def setUp(self):
        self.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx',
                                        'tests/fixtures/seq.fna', strict=False)
        self.model = wc_model_gen.ModelGenerator(
            knowledge_base=self.kb,
            component_generators=[prokaryote.ParametersGenerator]).run()

        self.model.id = 'prokaryote_compartment_tests'
        self._results_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._results_dir)

    def test_parameters(self):
        self.assertIsInstance(self.model.parameters.get_one(
            id='cellCycleLength'), wc_lang.core.Parameter)
        self.assertIsInstance(self.model.parameters.get_one(
            id='fractionDryWeight'), wc_lang.core.Parameter)
