""" Tests of RNA degradation submodel generation

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-07-24
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
                'ProteinDegradationSubmodelGenerator': {'beta': 1.}}}).run()

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

    def test_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='protein_degradation')
        
        modifier_species = model.observables.get_one(id='degrade_protease_obs').expression.species
        
        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)
            self.assertEqual(len(rl.expression.species), 3)
            self.assertEqual(set(rl.expression.species), set([
                i for i in rxn.get_reactants() if i not in modifier_species or i.species_type.id in rxn.id]))

        test_reaction = submodel.reactions.get_one(id='degradation_prot_gene_1_34')
        self.assertEqual(test_reaction.rate_laws[0].expression.expression, 
            'k_cat_degradation_prot_gene_1_34 * degrade_protease_obs * '
            '(prot_gene_1_34[c] / (prot_gene_1_34[c] + K_m_degradation_prot_gene_1_34_prot_gene_1_34 * Avogadro * volume_c)) * '
            '(atp[c] / (atp[c] + K_m_degradation_prot_gene_1_34_atp * Avogadro * volume_c)) * '
            '(h2o[c] / (h2o[c] + K_m_degradation_prot_gene_1_34_h2o * Avogadro * volume_c))')    

    def test_calibrate_submodel(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='protein_degradation')

        cytosol = kb.cell.compartments.get_one(id='c')
        test_species_type = kb.cell.species_types.get_one(id='prot_gene_1_1')
        test_species = test_species_type.species.get_one(compartment=cytosol)
        half_life = test_species_type.half_life
        mean_doubling_time = kb.cell.properties.get_one(id='mean_doubling_time').value
        protein_reactant_mean_concentration = kb.cell.concentrations.get_one(species=test_species).value
        degrase_mean_concentration = kb.cell.concentrations.get_one(
            species=kb.cell.species_types.get_one(id='prot_gene_1_34').species.get_one(compartment=cytosol)).value
        
        check_kcat = math.log(2) / half_life * protein_reactant_mean_concentration \
                    / degrase_mean_concentration / (0.5**3)

        self.assertAlmostEqual(model.parameters.get_one(id='k_cat_degradation_prot_gene_1_1').value, check_kcat, places=16)
