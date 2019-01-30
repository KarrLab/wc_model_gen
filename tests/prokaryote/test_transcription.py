""" Tests of transcription submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
from wc_utils.util.ontology import wcm_ontology
import math
import unittest
import wc_lang
import wc_kb
import wc_kb_gen


class TranscriptionSubmodelGeneratorTestCase(unittest.TestCase):

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
                                  prokaryote.TranscriptionSubmodelGenerator],
            options={'component': {
                'TranscriptionSubmodelGenerator': {
                    'rate_dynamics': 'phenomenological'}}}).run()

        cls.model_mechanistic = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=cls.kb,
            component_generators=[prokaryote.InitalizeModel,
                                  prokaryote.TranscriptionSubmodelGenerator],
            options={'component': {
                'TranscriptionSubmodelGenerator': {
                    'rate_dynamics': 'mechanistic'}}}).run()

    @classmethod
    def tearDownClass(cls):
        pass

    # Junk tests
    def test_submodels(self):
        model = self.model
        kb = self.kb

        submodel = model.submodels.get_one(id='transcription')
        self.assertIsInstance(submodel, wc_lang.Submodel)
        self.assertEqual(len(model.submodels), 2)

    def test_species(self):
        model = self.model
        kb = self.kb
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='transcription')

        for species in kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType):
            model_species = model.species_types.get_one(id=species.id)
            model_species_cytosol = model_species.species.get_one(compartment=cytosol)
            self.assertIsInstance(model_species, wc_lang.SpeciesType)
            self.assertIsInstance(model_species_cytosol, wc_lang.Species)

    def test_reactions(self):
        model = self.model
        kb = self.kb
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='transcription')

        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        # Check that number of RNAs = number of transcription reactions
        self.assertEqual(
            len(kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)),
            len(submodel.reactions))

        # Check that each reaction has the right number of participants
        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.participants), 10)

        # Check coeffs of reaction participants
        rnas = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rxn, rna in zip(submodel.reactions, rnas):

            self.assertEqual(
                + rxn.participants.get_one(species=atp).coefficient
                + rxn.participants.get_one(species=ctp).coefficient
                + rxn.participants.get_one(species=gtp).coefficient
                + rxn.participants.get_one(species=utp).coefficient,
                -rna.get_len())
            self.assertEqual(
                + rxn.participants.get_one(species=ppi).coefficient,
                rna.get_len())
            self.assertEqual(
                + rxn.participants.get_one(species=h2o).coefficient,
                -1)
            self.assertEqual(
                + rxn.participants.get_one(species=h).coefficient,
                rna.get_len() + 1)

    def test_phenom_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='transcription')

        for rxn in submodel.reactions:
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)
            self.assertEqual(len(rl.expression.species), 1)
            self.assertEqual(rl.expression.species[0].species_type.type, wcm_ontology['WCM:RNA']) # RNA
            self.assertIn(rl.expression.species[0], rxn.get_products())
