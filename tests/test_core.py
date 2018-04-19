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
        self.knowledge_base = wc_kb.KnowledgeBase()

    def test_ModelGenerator_constructor(self):
        generator = wc_model_gen.ModelGenerator(self.knowledge_base)

        self.assertEqual(generator.knowledge_base, self.knowledge_base)
        self.assertEqual(generator.components, ())
        self.assertEqual(generator.options, {})

    def test_ModelGenerator_run(self):
        generator = wc_model_gen.ModelGenerator(self.knowledge_base, options={
            'id': 'test_model',
        })
        model = generator.run()

        self.assertEqual(model.id, 'test_model')
        self.assertEqual(model.version, None)

    def test_ModelGenerator_run_with_components(self):
        class TestComponentGenerator1(wc_model_gen.ModelComponentGenerator):
            def run(self):
                self.model.compartments.create(id=self.options['compartment_id'])

        class TestComponentGenerator2(wc_model_gen.ModelComponentGenerator):
            def run(self): 
                self.model.compartments.create(name=self.options['compartment_name'])

        components = (
            TestComponentGenerator1, 
            TestComponentGenerator2,
            )
        options = {
            'id': 'test_model',
            'version': '0.0.1',
            'component': {
                'TestComponentGenerator1': {
                    'compartment_id': 'c'
                },
                'TestComponentGenerator2': {
                    'compartment_name': 'm'
                },
            },
        }
        gen = wc_model_gen.ModelGenerator(self.knowledge_base, components=components, options=options)
        model = gen.run()

        self.assertEqual(model.id, 'test_model')
        self.assertIsInstance(model.compartments.get_one(id='c'), wc_lang.Compartment)
        self.assertIsInstance(model.compartments.get_one(name='m'), wc_lang.Compartment)

    def test_ModelComponentGenerator_constructor(self):
        generator = wc_model_gen.ModelGenerator(self.knowledge_base)
        model = generator.run()

        class component_z_generator(wc_model_gen.ModelComponentGenerator):
            def run(self): pass

        component_z = component_z_generator(self.knowledge_base, model)

        self.assertEqual(component_z.knowledge_base, self.knowledge_base)
        self.assertEqual(component_z.model, model)
