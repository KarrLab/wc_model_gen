""" Tests of model generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""
from test.support import EnvironmentVarGuard
from wc_kb_gen import random
from wc_model_gen import prokaryote
import obj_model
import unittest
import wc_kb
import wc_lang
import tempfile
import wc_utils.util.string


class ModelGeneratorTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    @unittest.skip('fix charge/element imbalance')
    def test_submodels(self):
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_kb__DOT__io__DOT__strict', '0')
        with env:
            cls.kb = wc_kb.io.Reader().run('tests/fixtures/test_broken_kb.xlsx',
                                           'tests/fixtures/test_broken_seq.fna',
                                           )[wc_kb.KnowledgeBase][0]

        cls.model = prokaryote.ProkaryoteModelGenerator(knowledge_base=cls.kb).run()
        cytosol = self.model.compartments.get_one(id='c')
        extracellular_space = self.model.compartments.get_one(id='e')
        mean_doubling_time = self.model.parameters.get_one(id='mean_doubling_time')

        errors = obj_model.Validator().run(self.model, get_related=True)
        self.assertEqual(errors, None, msg=wc_utils.util.string.indent_forest(errors))
        self.assertEqual(5, len(self.model.submodels))
