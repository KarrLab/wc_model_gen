""" Tests of RNA degradation submodel generation

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-24
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
import math
import unittest
import wc_kb
import wc_kb_gen
import wc_lang


class ProteinDegradationSubmodelGeneratorTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_kb__DOT__io__DOT__strict', '0')
        with env:
            cls.kb = wc_kb.io.Reader().run('tests/fixtures/test_broken_kb.xlsx',
                                           'tests/fixtures/test_broken_seq.fna',
                                           )[wc_kb.KnowledgeBase][0]

        cls.model = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=cls.kb,
            component_generators=[prokaryote.InitalizeModel,
                                  prokaryote.ProteinDegradationSubmodelGenerator],
            options={'component': {
                'ProteinDegradationSubmodelGenerator': {
                    'rate_dynamics': 'phenomenological'}}}).run()

        cls.model_mechanistic = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=cls.kb,
            component_generators=[prokaryote.InitalizeModel,
                                  prokaryote.ProteinDegradationSubmodelGenerator],
            options={'component': {
                'ProteinDegradationSubmodelGenerator': {
                    'rate_dynamics': 'mechanistic'}}}).run()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_species(self):
        model = self.model
        kb = self.kb
        cell = self.kb.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='protein_degradation')

        # check reactions generated
        self.assertEqual(len(submodel.reactions),
                         len(cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)))

        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        adp = model.species_types.get_one(id='adp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)

        # check species types and species generated
        for species in kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType):
            model_species = model.species_types.get_one(id=species.id)
            model_species_cytosol = model_species.species.get_one(compartment=cytosol)
            self.assertIsInstance(model_species, wc_lang.SpeciesType)
            self.assertIsInstance(model_species_cytosol, wc_lang.Species)

        for rxn in submodel.reactions:
            self.assertEqual(rxn.participants.get_one(species=atp).coefficient, -1)
            self.assertEqual(rxn.participants.get_one(species=adp).coefficient, 1)
            self.assertEqual(rxn.participants.get_one(species=pi).coefficient, 1)

        # TODO: add coutn of AAs
        # self.assertEqual(submodel.reactions[0].participants.get_one(
        #    species=aa_species).coefficient, prots[0].get_seq().count('C'))

    def test_phenom_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='protein_degradation')

        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)
            self.assertEqual(len(rl.expression.species), 1)
            self.assertEqual(rl.expression.species[0].species_type.type, wc_lang.SpeciesTypeType.protein)
            self.assertIn(rl.expression.species[0], rxn.get_reactants())

    def test_mechanistic_rate_laws(self):
        model = self.model_mechanistic
        kb = self.kb
        submodel = model.submodels.get_one(id='protein_degradation')

        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)

            # TODO:
            self.assertEqual(rxn.get_modifiers(), [])

            k_cat_value = rl.expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value
            self.assertIsInstance(k_cat_value, float)
            self.assertFalse(math.isnan(k_cat_value))

            # Check that reactants participate in rate law
            self.assertEqual(set(rxn.get_reactants()).difference(set(rl.expression.species)), set())
