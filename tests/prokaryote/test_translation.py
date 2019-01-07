""" Testing Translation Submodel Generator

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-23
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
import math
import unittest
import wc_kb
import wc_lang


class TranslationSubmodelGeneratorTestCase(unittest.TestCase):

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
                                  prokaryote.TranslationSubmodelGenerator],
            options={'component': {
                'TranslationSubmodelGenerator': {
                    'rate_dynamics': 'phenomenological'}}}).run()

        cls.model_mechanistic = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=cls.kb,
            component_generators=[prokaryote.InitalizeModel,
                                  prokaryote.TranslationSubmodelGenerator],
            options={'component': {
                'TranslationSubmodelGenerator': {
                    'rate_dynamics': 'mechanistic'}}}).run()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_submodels(self):
        kb = self.kb
        model = self.model
        model_mechanistic = self.model_mechanistic

        submodel = model.submodels.get_one(id='translation')
        self.assertIsInstance(submodel, wc_lang.Submodel)
        self.assertEqual(len(model.submodels), 2)

        submodel = model_mechanistic.submodels.get_one(id='translation')
        self.assertIsInstance(submodel, wc_lang.Submodel)
        self.assertEqual(len(model_mechanistic.submodels), 2)

    def test_species(self):
        model = self.model
        kb = self.kb
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='translation')

        for species in kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType):
            model_species = model.species_types.get_one(id=species.id)
            model_species_cytosol = model_species.species.get_one(compartment=cytosol)
            self.assertIsInstance(model_species, wc_lang.SpeciesType)
            self.assertIsInstance(model_species_cytosol, wc_lang.Species)

    def test_reactions(self):
        model = self.model
        kb = self.kb
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='translation')

        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        gdp = model.species_types.get_one(id='gdp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)
        ribosome = model.observables.get_one(id='ribosome_obs').expression.species[0]
        initiation_factors = model.observables.get_one(id='translation_init_factors_obs').expression.species[0]
        elongation_factors = model.observables.get_one(id='translation_elongation_factors_obs').expression.species[0]
        release_factors = model.observables.get_one(id='translation_release_factors_obs').expression.species[0]

        # Check that number of RNAs = number of transcription reactions
        self.assertEqual(
            len(kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)),
            len(submodel.reactions))

        # Check that each reaction has the min or more number of participants
        for rxn in submodel.reactions:
            self.assertTrue(len(rxn.participants) > 5)

        # Check coeffs of reaction participants
        # TODO: add assertions about the number of participating tRNAs
        prots_kb = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for rxn, prot_kb in zip(submodel.reactions, prots_kb):
            prot_model = model.species_types.get_one(id=prot_kb.id)
            length = len(prot_kb.get_seq())

            self.assertEqual(rxn.participants.get_one(species=gtp).coefficient, -(length+2))
            self.assertEqual(rxn.participants.get_one(species=gdp).coefficient, (length+2))
            self.assertEqual(rxn.participants.get_one(species=pi).coefficient, 2*length)

            """
            Need to customize assertions for translation reactions that produce:
            ribosome
            initiation_factors
            elongation_factors
            release_factors

            self.assertEqual(len(rxn.participants.get(species=ribosome)), 2)
            self.assertEqual(abs(rxn.participants.get(species=ribosome)[0].coefficient), 1)
            self.assertEqual(abs(rxn.participants.get(species=ribosome)[1].coefficient), 1)

            self.assertEqual(len(rxn.participants.get(species=initiation_factors)), 2)
            self.assertEqual(abs(rxn.participants.get(species=initiation_factors)[0].coefficient), 1)
            self.assertEqual(abs(rxn.participants.get(species=initiation_factors)[1].coefficient), 1)

            self.assertEqual(len(rxn.participants.get(species=initiation_factors))=2)
            self.assertEqual(abs(rxn.participants.get(species=initiation_factors)[0].coefficient)=1)
            self.assertEqual(abs(rxn.participants.get(species=initiation_factors)[1].coefficient)=1)

            self.assertEqual(rxn.participants.get_one(species=initiation_factors).coefficient, -(length+2))
            self.assertEqual(rxn.participants.get_one(species=elongation_factors).coefficient, -length)
            self.assertEqual(rxn.participants.get_one(species=release_factors).coefficient, (length+2))
            """

    def test_phenom_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='translation')

        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)
            self.assertEqual(len(rl.expression.species), 1)
            self.assertEqual(rl.expression.species[0].species_type.type, wc_lang.SpeciesTypeType.protein)
            self.assertIn(rl.expression.species[0], rxn.get_products())

    def test_mechanistic_rate_laws(self):
        model = self.model_mechanistic
        kb = self.kb
        submodel = model.submodels.get_one(id='translation')

        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, 1)
            self.assertEqual(rxn.get_modifiers(), [])

            k_cat_value = rl.expression.parameters.get_one(type=wc_lang.ParameterType.k_cat).value
            self.assertIsInstance(k_cat_value, float)
            self.assertFalse(math.isnan(k_cat_value))

            # Check that reactants participate in rate law
            self.assertEqual(set(rxn.get_reactants()).difference(set(rl.expression.species)), set())
