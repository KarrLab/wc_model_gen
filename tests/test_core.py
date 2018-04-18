""" Tests of the model generator

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-02-07
:Copyright: 2018, Karr Lab
:License: MIT
"""

import model_generator
import unittest
import wc_kb

class TestModelGenerator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.knowledge_base = kb = wc_kb.KnowledgeBase()

    def test_ModelGenerator_constructor(self):
        generator = model_generator.ModelGenerator(self.knowledge_base)

        self.assertEqual(generator.knowledge_base, self.knowledge_base)
        self.assertEqual(generator.version, '0.0.1')
        self.assertEqual(generator.components, ())
        self.assertEqual(generator.options, {})

    def test_ModelGenerator_run(self):
        generator = model_generator.ModelGenerator(self.knowledge_base)
        model = generator.run()

        self.assertEqual(model.id, 'test_model')

    def test_ModelComponentGenerator_constructor(self):
        generator = model_generator.ModelGenerator(self.knowledge_base)
        model = generator.run()

        class component_z_generator(model_generator.ModelComponentGenerator):
            def run(self): pass

        component_z = component_z_generator(self.knowledge_base, model)

        self.assertEqual(component_z.knowledge_base, self.knowledge_base)
        self.assertEqual(component_z.model, model)
