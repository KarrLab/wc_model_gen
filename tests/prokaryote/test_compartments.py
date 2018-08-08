import shutil
import os
import unittest
import tempfile
import wc_kb
import wc_lang
import wc_model_gen
import wc_model_gen.prokaryote as prokaryote


class CompartmentGeneratorTestCase(unittest.TestCase):

    def setUp(self):
        self.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx',
                                        'tests/fixtures/seq.fna', strict=False)
        self.model = wc_model_gen.ModelGenerator(
            knowledge_base=self.kb,
            component_generators=[prokaryote.CompartmentsGenerator,
                                  prokaryote.ParametersGenerator,
                                  prokaryote.MetabolismSubmodelGenerator]).run()

        self.model.id = 'prokaryote_compartment_tests'
        self._results_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._results_dir)

    def test_compartments(self):
        cytosol = self.model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')
        ec = self.model.compartments.get_one(id='e')
        self.assertEqual(ec.name, 'extracellular space')
