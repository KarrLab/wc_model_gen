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

    def test_ModelGenerator_run_with_submodels(self):

        class TestComponentGenerator1(wc_model_gen.ModelComponentGenerator):
            def run(self):
                self.model.compartments.create(id=self.options['compartment_id'], name='Cytosol')

        class TestSubmodelGenerator2(wc_model_gen.SubmodelGenerator):
            def gen_species(self):
                self.model.species_types.create(id=self.options['species_type_id'], name='ATP')

            def gen_phenomenological_rates(self):
                pass

        class TestSubmodelGenerator3(wc_model_gen.SubmodelGenerator):
            def gen_species(self):
                self.model.species_types.create(id='gtp', name=self.options['species_type_name'])

            def gen_phenomenological_rates(self):
                pass

        component_generators = [TestComponentGenerator1,
                                TestSubmodelGenerator2,
                                TestSubmodelGenerator3]

        options = {
            'id': 'test_model',
            'version': '0.0.1',
            'component': {
                'TestComponentGenerator1': {
                    'compartment_id': 'c'
                },
                'TestSubmodelGenerator2': {
                    'species_type_id': 'atp'
                },
                'TestSubmodelGenerator3': {
                    'species_type_name': 'GTP'
                },
            },
        }

        kb = wc_kb.core.KnowledgeBase()
        model = wc_model_gen.ModelGenerator(kb, component_generators=component_generators, options=options).run()

        self.assertIsInstance(model.compartments.get_one(id='c'), wc_lang.core.Compartment)
        self.assertIsInstance(model.species_types.get_one(id='atp'), wc_lang.core.SpeciesType)
        self.assertIsInstance(model.species_types.get_one(name='GTP'), wc_lang.core.SpeciesType)

    def test_SubmodelGenerator(self):
        model = wc_model_gen.ModelGenerator(self.knowledge_base).run()

        class component_z_generator(wc_model_gen.SubmodelGenerator):
            #def gen_compartments(self): pass

            def gen_species(self): pass

            def gen_reactions(self): pass

            def gen_rate_laws(self): pass

        component_z = component_z_generator(self.knowledge_base, model)

        self.assertEqual(component_z.knowledge_base, self.knowledge_base)
        self.assertEqual(component_z.model, model)
