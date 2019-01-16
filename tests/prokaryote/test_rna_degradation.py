""" Tests of RNA degradation submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
from wc_utils.util.ontology import wcm_ontology
import math
import unittest
import wc_kb
import wc_kb_gen
import wc_lang


class RnaDegradationSubmodelGeneratorTestCase(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_kb__DOT__io__DOT__strict', '0')
        with env:
            cls.kb = wc_kb.io.Reader().run('tests/fixtures/test_broken_kb.xlsx',
                                           'tests/fixtures/test_broken_seq.fna',
                                            )[wc_kb.KnowledgeBase][0]

        cls.model = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.RnaDegradationSubmodelGenerator],
                        options = {'component': {
                             'RnaDegradationSubmodelGenerator': {
                               'rate_dynamics': 'phenomenological'}}}).run()

        cls.model_mechanistic = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.RnaDegradationSubmodelGenerator],
                        options = {'component': {
                             'RnaDegradationSubmodelGenerator': {
                               'rate_dynamics': 'mechanistic'}}}).run()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_species(self):
        model = self.model
        kb = self.kb
        cell = self.kb.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='rna_degradation')

        # check reactions generated
        self.assertEqual(len(submodel.reactions),
                         len(cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)))

        # check species types and species generated
        for species in kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType):
            model_species = model.species_types.get_one(id=species.id)
            model_species_cytosol = model_species.species.get_one(compartment=cytosol)
            self.assertIsInstance(model_species, wc_lang.SpeciesType)
            self.assertIsInstance(model_species_cytosol, wc_lang.Species)


        amp = model.species_types.get_one(id='amp').species.get_one(compartment=cytosol)
        cmp = model.species_types.get_one(id='cmp').species.get_one(compartment=cytosol)
        gmp = model.species_types.get_one(id='gmp').species.get_one(compartment=cytosol)
        ump = model.species_types.get_one(id='ump').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        # Check coeffs of reaction participants
        rnas = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rxn, rna in zip(submodel.reactions, rnas):
            self.assertEqual(
                + rxn.participants.get_one(species=amp).coefficient
                + rxn.participants.get_one(species=cmp).coefficient
                + rxn.participants.get_one(species=gmp).coefficient
                + rxn.participants.get_one(species=ump).coefficient,
                rna.get_len())
            self.assertEqual(
                + rxn.participants.get_one(species=h2o).coefficient,
                -(rna.get_len() - 1))
            self.assertEqual(
                + rxn.participants.get_one(species=h).coefficient,
                rna.get_len() - 1)

    def test_phenom_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='rna_degradation')

        for rxn in submodel.reactions:

            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)
            self.assertEqual(len(rl.expression.species), 1)
            self.assertEqual(rl.expression.species[0].species_type.type, wcm_ontology['WCM:RNA']) # RNA
            self.assertIn(rl.expression.species[0], rxn.get_reactants())

    def test_mechanistic_rate_laws(self):
        model = self.model_mechanistic
        kb = self.kb
        submodel = model.submodels.get_one(id='rna_degradation')

        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)
            
            self.assertEqual(rxn.get_modifiers(), [])

            k_cat_value = rl.expression.parameters.get_one(type=wcm_ontology['WCM:k_cat']).value
            self.assertIsInstance(k_cat_value, float)
            self.assertFalse(math.isnan(k_cat_value))

            # Check that reactants participate in rate law
            self.assertEqual(set(rxn.get_reactants()).difference(set(rl.expression.species)), set())
