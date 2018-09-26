""" Testing Translation Submodel Generator

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-23
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_lang
import wc_kb


class TranslationSubmodelGeneratorTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/min_kb.xlsx',
                                       'tests/fixtures/min_kb_seq.fna',
                                        strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.TranslationSubmodelGenerator],
                        options = {'component': {
                             'TranslationSubmodelGenerator': {
                               'rate_dynamics': 'phenomenological'}}}).run()

        cls.model_mechanistic = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.TranslationSubmodelGenerator],
                        options = {'component': {
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
        self.assertIsInstance(submodel, wc_lang.core.Submodel)
        self.assertEqual(len(model.submodels), 1)

        submodel = model_mechanistic.submodels.get_one(id='translation')
        self.assertIsInstance(submodel, wc_lang.core.Submodel)
        self.assertEqual(len(model_mechanistic.submodels), 1)

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

        #Check that number of RNAs = number of transcription reactions
        self.assertEqual(
            len(kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)),
            len(submodel.reactions))

        # Check coeffs of reaction participants
        """
        prots_kb = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for rxn, prot_kb in zip(submodel.reactions, prots_kb):

            prot_model = model.species_types.get_one(id=prot_kb)
            length = len(prot_kb.get_seq())

            self.assertEqual(rxn.participants.get_one(species=prot_model.species[0]).coefficient, 1)
            self.assertEqual(rxn.participants.get_one(species=gtp).coefficient, -(length+2))
            self.assertEqual(rxn.participants.get_one(species=gdp).coefficient, (length+2))
            self.assertEqual(rxn.participants.get_one(species=pi).coefficient, 2*length)
        """

    @unittest.skip
    def test_rate_laws(self):
        pass #TODO
