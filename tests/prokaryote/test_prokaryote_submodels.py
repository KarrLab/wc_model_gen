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
        shutil.rmtree(cls._results_dir)

    def test_submodel_added(self):
        submodel_ids = []
        for submodel in self.model.get_submodels():
            submodel_ids.append(submodel.id)

        self.assertTrue('transcription' in submodel_ids)

    # Check if transcription proceeds without triphosphate molecules
    def test_are_triphosphates_bottleneck_to_transcription(self):
        pass
