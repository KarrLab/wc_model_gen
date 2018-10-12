""" Tests of RNA degradation submodel generation

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-24
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb_gen
import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_lang
import wc_kb
import math


class ProteinDegradationSubmodelGeneratorTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/test_broken_kb.xlsx',
                                       'tests/fixtures/test_broken_seq.fna',
                                        strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.ProteinDegradationSubmodelGenerator],
                        options = {'component': {
                             'ProteinDegradationSubmodelGenerator': {
                               'rate_dynamics': 'phenomenological'}}}).run()

        cls.model_mechanistic = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.ProteinDegradationSubmodelGenerator],
                        options = {'component': {
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
        #self.assertEqual(submodel.reactions[0].participants.get_one(
        #    species=aa_species).coefficient, prots[0].get_seq().count('C'))

    def test_phenom_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='protein_degradation')

        for rxn in submodel.reactions:

            self.assertEqual(len(rxn.rate_laws), 1)
            self.assertIsInstance(rxn.rate_laws[0], wc_lang.core.RateLaw)
            self.assertEqual(rxn.rate_laws[0].direction, 1)
            self.assertEqual(len(rxn.rate_laws[0].equation.modifiers), 1)

            # Check that RNA produced is modifier
            match = 0
            for participant in rxn.participants:
                if participant.species == rxn.rate_laws[0].equation.modifiers[0]:
                    match = 1
                    break

            self.assertEqual(match, 1)

    def test_mechanistic_rate_laws(self):
        model = self.model_mechanistic
        kb = self.kb
        submodel = model.submodels.get_one(id='protein_degradation')

        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            self.assertIsInstance(rxn.rate_laws[0], wc_lang.core.RateLaw)
            self.assertEqual(rxn.rate_laws[0].direction, 1)
            #print(len(rxn.rate_laws[0].equation.modifiers))
            #TODO:
            self.assertTrue(len(rxn.rate_laws[0].equation.modifiers) >= 3)

            self.assertIsInstance(rxn.rate_laws[0].k_cat, float)
            self.assertFalse(math.isnan(rxn.rate_laws[0].k_cat))

            # Check that participants are modifiers
            for participant in rxn.participants:
                if participant.coefficient < 0:
                    self.assertTrue(participant.species in rxn.rate_laws[0].equation.modifiers)
