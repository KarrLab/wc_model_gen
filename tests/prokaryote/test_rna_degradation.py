""" Tests of RNA degradation submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
from wc_onto import onto as wc_ontology
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
            cls.kb = wc_kb.io.Reader().run('tests/fixtures/min_model_kb.xlsx',
                                           'tests/fixtures/min_model_seq.fna',
                                            )[wc_kb.KnowledgeBase][0]

        cls.model = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.RnaDegradationSubmodelGenerator],
                        options = {'component': {
                             'RnaDegradationSubmodelGenerator': {'beta': 1.}}}).run()

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
                         len(cell.species_types.get(__type=wc_kb.prokaryote.RnaSpeciesType)))

        # check species types and species generated
        for species in kb.cell.species_types.get(__type=wc_kb.prokaryote.RnaSpeciesType):
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
        rnas = kb.cell.species_types.get(__type=wc_kb.prokaryote.RnaSpeciesType)
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

    def test_rate_laws(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='rna_degradation')

        modifier_species = model.observables.get_one(id='degrade_rnase_obs').expression.species

        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertIsInstance(rl, wc_lang.RateLaw)
            self.assertEqual(rl.direction, wc_lang.RateLawDirection.forward)
            self.assertEqual(len(rl.expression.species), 2)
            self.assertEqual(set(rl.expression.species), set([i for i in rxn.get_reactants() if i not in modifier_species]))

        test_reaction = submodel.reactions.get_one(id='degradation_rna_tu_1_1')
        self.assertEqual(test_reaction.rate_laws[0].expression.expression,
            'k_cat_degradation_rna_tu_1_1 * degrade_rnase_obs * '
            '(rna_tu_1_1[c] / (rna_tu_1_1[c] + K_m_degradation_rna_tu_1_1_rna_tu_1_1 * Avogadro * volume_c)) * '
            '(h2o[c] / (h2o[c] + K_m_degradation_rna_tu_1_1_h2o * Avogadro * volume_c))')

    def test_calibrate_submodel(self):
        model = self.model
        kb = self.kb
        submodel = model.submodels.get_one(id='rna_degradation')

        cytosol = kb.cell.compartments.get_one(id='c')
        test_species_type = kb.cell.species_types.get_one(id='rna_tu_1_1')
        test_species = test_species_type.species.get_one(compartment=cytosol)
        half_life = test_species_type.properties.get_one(property='half_life').get_value()
        mean_doubling_time = kb.cell.parameters.get_one(id='mean_doubling_time').value
        rna_reactant_mean_concentration = kb.cell.concentrations.get_one(species=test_species).value
        degrase_mean_concentration = kb.cell.concentrations.get_one(
            species=kb.cell.species_types.get_one(id='prot_gene_1_11').species.get_one(compartment=cytosol)).value

        check_kcat = math.log(2) / half_life * rna_reactant_mean_concentration \
                    / degrase_mean_concentration / (0.5**2)

        self.assertAlmostEqual(model.parameters.get_one(id='k_cat_degradation_rna_tu_1_1').value, check_kcat, places=16)
