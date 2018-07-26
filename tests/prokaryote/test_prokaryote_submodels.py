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

class ProkaryoteTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx','tests/fixtures/seq.fna', strict = False)
        cls.model = wc_model_gen.ModelGenerator(
                            knowledge_base=cls.kb,
                            component_generators=[prokaryote.CompartmentsGenerator,
                                                  prokaryote.ParametersGenerator,
                                                  prokaryote.MetaboliteSpeciesGenerator]).run()

        cls.model.id = 'prokaryote_submodels_tests'
        cls._results_dir = tempfile.mkdtemp()
        cls.model.version = '0.0.1'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._results_dir)

class TranscriptionSubmodelTestCase(ProkaryoteTestCase):

    @classmethod
    def setUpClass(cls):
        super(TranscriptionSubmodelTestCase, cls).setUpClass()
        prokaryote.TranscriptionSubmodelGenerator(cls.kb, cls.model).run()

    def test_submodel_added(self):
        submodel_ids = []
        for submodel in self.model.get_submodels():
            submodel_ids.append(submodel.id)

        self.assertTrue('transcription' in submodel_ids)

class TranslationSubmodelTestCase(ProkaryoteTestCase):

    @classmethod
    def setUpClass(cls):
        super(TranslationSubmodelTestCase, cls).setUpClass()
        prokaryote.TranslationSubmodelGenerator(cls.kb, cls.model).run()

    def test_submodel_added(self):
        submodel_ids = []
        for submodel in self.model.get_submodels():
            submodel_ids.append(submodel.id)

        self.assertTrue('translation' in submodel_ids)

class DegradationSubmodelTestCase(ProkaryoteTestCase):

    @classmethod
    def setUpClass(cls):
        super(DegradationSubmodelTestCase, cls).setUpClass()
        prokaryote.TranscriptionSubmodelGenerator(cls.kb, cls.model).run()

    def test_submodel_added(self):
        submodel_ids = []
        for submodel in self.model.get_submodels():
            submodel_ids.append(submodel.id)

        self.assertTrue('degradation' in submodel_ids)
