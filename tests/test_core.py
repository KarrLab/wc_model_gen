""" Tests of the model generator

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen
import unittest
import wc_kb
import wc_lang


class TestModelGenerator(unittest.TestCase):
    def setUp(self):
        self.knowledge_base = wc_kb.core.KnowledgeBase()

    def test_ModelGenerator_init(self):
        generator = wc_model_gen.ModelGenerator(self.knowledge_base)

        self.assertEqual(generator.knowledge_base, self.knowledge_base)
        self.assertEqual(generator.component_generators, [])
        self.assertEqual(generator.options, {'id': None, 'name': None, 'version': None})

    def test_ModelGenerator_run(self):
        generator = wc_model_gen.ModelGenerator(self.knowledge_base, options={
            'id': 'test_model'})

        model = generator.run()

        self.assertEqual(model.id, 'test_model')
        self.assertEqual(model.version, None)

        generator = wc_model_gen.ModelGenerator(self.knowledge_base, options={
            'id': 'test_model',
            'version': '0.0.1',
        })
        model = generator.run()

        self.assertEqual(model.id, 'test_model')
        self.assertEqual(model.version, '0.0.1')

    def test_SubmodelGenerator(self):
        model = wc_model_gen.ModelGenerator(self.knowledge_base).run()

        class component_z_generator(wc_model_gen.SubmodelGenerator):
            def gen_reactions(self): pass

            def gen_rate_laws(self): pass

        component_z = component_z_generator(self.knowledge_base, model)

        self.assertEqual(component_z.knowledge_base, self.knowledge_base)
        self.assertEqual(component_z.model, model)
