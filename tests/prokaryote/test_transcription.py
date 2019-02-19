""" Tests of transcription submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
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
                'TranscriptionSubmodelGenerator': {'beta': 1.}}}).run()

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
            self.assertEqual(len(rxn.participants), 6)

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
                rna.get_len() -1)

    def test_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='transcription')
        
        for rxn in submodel.reactions:
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)
            self.assertEqual(len(rl.expression.species), 4)            
            self.assertEqual(set(rl.expression.species), set(rxn.get_reactants()))

        test_reaction = submodel.reactions.get_one(id='transcription_rna_tu_1_1')
        self.assertEqual(test_reaction.rate_laws[0].expression.expression, 
            'k_cat_transcription_rna_tu_1_1 * rna_polymerase_obs * '
            '(atp[c] / (atp[c] + K_m_transcription_rna_tu_1_1_atp * Avogadro * volume_c)) * '
            '(ctp[c] / (ctp[c] + K_m_transcription_rna_tu_1_1_ctp * Avogadro * volume_c)) * '
            '(gtp[c] / (gtp[c] + K_m_transcription_rna_tu_1_1_gtp * Avogadro * volume_c)) * '
            '(utp[c] / (utp[c] + K_m_transcription_rna_tu_1_1_utp * Avogadro * volume_c))')    

    def test_calibrate_submodel(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='transcription')
        
        cytosol = kb.cell.compartments.get_one(id='c')
        test_species_type = kb.cell.species_types.get_one(id='rna_tu_1_1')
        test_species = test_species_type.species.get_one(compartment=cytosol)
        half_life = test_species_type.half_life
        mean_doubling_time = kb.cell.properties.get_one(id='mean_doubling_time').value
        rna_mean_concentration = kb.cell.concentrations.get_one(species=test_species).value
        prot_mean_concentration = kb.cell.concentrations.get_one(
            species=kb.cell.species_types.get_one(id='prot_gene_1_31').species.get_one(compartment=cytosol)).value
        
        check_kcat = math.log(2) * (1. / half_life + 1. / mean_doubling_time) * rna_mean_concentration \
                    / prot_mean_concentration / (0.5**4)

        self.assertAlmostEqual(model.parameters.get_one(id='k_cat_transcription_rna_tu_1_1').value, check_kcat, places=16) 
