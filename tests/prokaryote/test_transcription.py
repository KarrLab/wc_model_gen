""" Tests of transcription submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_model_gen import prokaryote
import math
import unittest
import wc_kb
import wc_kb_gen
import wc_lang

class TranscriptionSubmodelGeneratorTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/test_broken_kb.xlsx',
                                       'tests/fixtures/test_broken_seq.fna',
                                        strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.TranscriptionSubmodelGenerator],
                        options = {'component': {
                             'TranscriptionSubmodelGenerator': {
                               'rate_dynamics': 'phenomenological'}}}).run()

        cls.model_mechanistic = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.TranscriptionSubmodelGenerator],
                        options = {'component': {
                             'TranscriptionSubmodelGenerator': {
                               'rate_dynamics': 'mechanistic'}}}).run()

    @classmethod
    def tearDownClass(cls):
        pass

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
        h   = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        #Check that number of RNAs = number of transcription reactions
        self.assertEqual(
            len(kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)),
            len(submodel.reactions))

        # Check that each reaction has the right number of participants
        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.participants),10)

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
                rna.get_len() - 1)
            self.assertEqual(
                + rxn.participants.get_one(species=h).coefficient,
                -(rna.get_len() - 1))

    def test_phenom_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='transcription')

        for rxn in submodel.reactions:

            self.assertEqual(len(rxn.rate_laws), 1)
            self.assertIsInstance(rxn.rate_laws[0], wc_lang.RateLaw)
            self.assertEqual(rxn.rate_laws[0].direction, 1)
            self.assertEqual(len(rxn.rate_laws[0].equation.modifiers), 1)

            # Check that RNA produced is modifier
            match = 0
            for participant in rxn.participants:
                if participant.species == rxn.rate_laws[0].equation.modifiers[0]:
                    match = 1
                    break

            self.assertEqual(match, 1)
