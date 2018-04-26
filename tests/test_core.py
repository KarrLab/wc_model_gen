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


class TestInitalizeModel(unittest.TestCase):
    def setUp(self):
        self.knowledge_base = wc_kb.KnowledgeBase()

    def test_InitalizeModel(self):
        generator = wc_model_gen.InitalizeModel(self.knowledge_base)

        self.assertEqual(generator.knowledge_base, self.knowledge_base)
        self.assertEqual(generator.components, [])
        self.assertEqual(generator.options, {})

    def test_InitalizeModel_run(self):
        generator = wc_model_gen.InitalizeModel(self.knowledge_base, options={
            'id': 'test_model',
        })
        model = generator.run()

        self.assertEqual(model.id, 'test_model')
        self.assertEqual(model.version, '0.0.1')

    def test_InitalizeModel_run_with_components(self):
        class TestComponentGenerator1(wc_model_gen.SubmodelGenerator):
            def run(self):
                self.model.compartments.create(id='c', name='Cytosol')
            def generate_species(self): pass
            def generate_reactions(self): pass
            def generate_rate_laws(self): pass


        class TestComponentGenerator2(wc_model_gen.SubmodelGenerator):
            def run(self):
                self.model.compartments.create(id='m', name='Membrane')
            def generate_species(self): pass
            def generate_reactions(self): pass
            def generate_rate_laws(self): pass

        components = [
            TestComponentGenerator1,
            TestComponentGenerator2,
            ]

        """ WIP
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
        """
        kb = wc_kb.KnowledgeBase()
        print(type(kb))
        model = wc_model_gen.InitalizeModel(kb,components=components).run()

        print(model.compartments[0].id)

        self.assertIsInstance(model.compartments.get_one(id='c'), wc_lang.core.Compartment)
        self.assertIsInstance(model.compartments.get_one(name='Membrane'), wc_lang.core.Compartment)

    def test_SubmodelGenerator(self):
        model = wc_model_gen.InitalizeModel(self.knowledge_base).run()

        class component_z_generator(wc_model_gen.SubmodelGenerator):
            def generate_species(self): pass
            def generate_reactions(self): pass
            def generate_rate_laws(self): pass

        component_z = component_z_generator(self.knowledge_base, model)

        self.assertEqual(component_z.knowledge_base, self.knowledge_base)
        self.assertEqual(component_z.model, model)
