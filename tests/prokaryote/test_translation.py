""" Testing Translation Submodel Generator

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-07-23
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
from wc_onto import onto as wc_ontology
import math
import scipy.constants
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
                'TranslationSubmodelGenerator': {'beta': 1.}}}).run()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_submodels(self):
        kb = self.kb
        model = self.model
        
        submodel = model.submodels.get_one(id='translation')
        self.assertIsInstance(submodel, wc_lang.Submodel)
        self.assertEqual(len(model.submodels), 2)

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

    def test_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='translation')

        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, 1)

        #TODO: Update test once translation module is complete
        test_reaction = submodel.reactions.get_one(id='translation_prot_gene_1_1')
        self.assertEqual(len(test_reaction.rate_laws[0].expression.parameters), 7) # kcat + Avogadro + 5 Kms (4 aa + 1 gtp)
        self.assertEqual(len(test_reaction.rate_laws[0].expression.observables), 8) # 3 factors + 4 tRNAs + ribosome
        self.assertEqual(len(test_reaction.rate_laws[0].expression.species), 5) # 4 aa + 1 gtp
        self.assertEqual(len(test_reaction.rate_laws[0].expression.functions), 1) # volume
                        
    def test_calibrate_submodel(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='transcription')
        
        #TODO: Update test once translation module is complete        
        cytosol = kb.cell.compartments.get_one(id='c')
        test_species_type = kb.cell.species_types.get_one(id='prot_gene_1_1')
        test_species = test_species_type.species.get_one(compartment=cytosol)
        half_life = test_species_type.half_life
        mean_doubling_time = kb.cell.parameters.get_one(id='mean_doubling_time').value

        Avogadro = scipy.constants.Avogadro
        volume = model.compartments.get_one(id='c').mean_init_volume
        
        protein_mean_concentration = kb.cell.concentrations.get_one(species=test_species).value * volume * Avogadro
        
        modifier_species_type = ['prot_gene_1_10', 'prot_gene_1_13', 'prot_gene_1_29', 
                                'rna_tu_1_33', 'rna_tu_1_28', 'rna_tu_1_24', 'rna_tu_1_7',
                                'ribosome'] # 3 factors + 4 tRNAs + ribosome
        observables_products = 1.
        for i in modifier_species_type:
            observables_products *= kb.cell.concentrations.get_one(
                species=kb.cell.species_types.get_one(id=i).species.get_one(compartment=cytosol)).value * volume * Avogadro
        
        check_kcat = math.log(2) * (1. / half_life + 1. / mean_doubling_time) * protein_mean_concentration \
                    / observables_products / (0.5**5) # **5 because 4 aa + 1 gtp

        self.assertAlmostEqual(model.parameters.get_one(id='k_cat_translation_prot_gene_1_1').value, check_kcat, places=30)
              