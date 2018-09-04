""" Tests of metabolism submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""
from wc_model_gen.prokaryote import metabolism
import unittest
import wc_lang
import wc_kb


class MetabolismSubmodelGeneratorTestCase(unittest.TestCase):

    def setUp(self):

        self.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx',
                                        'tests/fixtures/seq.fna', strict=False)
        self.model = wc_lang.Model()
        metabolism.MetabolismSubmodelGenerator(
            self.kb, self.model, options={}).run()
        self.submodel = self.model.submodels.get_one(id='metabolism')

    def testModelGen(self):
        self.assertIsInstance(self.submodel, wc_lang.Submodel)

    def testSpecies(self):

        # check parameters generated
        self.assertEqual(self.model.parameters.get_one(
            id='fraction_dry_weight').value, 0.7)

        # check species types and species generated
        cytosol = self.model.compartments.get(id='c')[0]

        for species in self.kb.cell.species_types.get(__type=wc_kb.MetaboliteSpeciesType):
            model_species = self.model.species_types.get_one(id=species.id)
            model_species_cytosol = model_species.species.get_one(
                compartment=cytosol)
            self.assertIsInstance(model_species, wc_lang.SpeciesType)
            self.assertIsInstance(model_species_cytosol, wc_lang.Species)
            self.assertEqual(
                model_species_cytosol.concentration.units, wc_lang.ConcentrationUnit.M)
            # self.assertEqual(
            #    model_species_cytosol.concentration, species.concentration)

    def testReactions(self):
        # aggregate properties
        kb_reactions = self.kb.cell.reactions
        lang_reactions = self.submodel.reactions
        self.assertEqual(len(kb_reactions), len(lang_reactions))
        self.assertEqual(set([r.id for r in kb_reactions]), set([r.id for r in lang_reactions]))

    def get_species(self, species_type_id, compartment_id):
        lang_species_type = self.model.species_types.get_one(
            id=species_type_id)
        lang_species = lang_species_type.species.get_one(
            compartment=self.model.compartments.get_one(id=compartment_id))
        return lang_species

    def test_o2_degradation_rxn(self):
        lang_O2_degradation_rxn = self.submodel.reactions.get_one(id='O2_degradation')
        self.assertEqual(lang_O2_degradation_rxn.name, 'O2 degradation')
        self.assertEqual(lang_O2_degradation_rxn.reversible, False)
        kb_O2_degradation = self.kb.cell.reactions.get_one(id='O2_degradation')
        expected_participants = set()
        for part in kb_O2_degradation.participants:
            lang_species = self.get_species(part.species.species_type.id, part.species.compartment.id)
            expected_participants.add(lang_species.species_coefficients.get_one(
                coefficient=part.coefficient))
        self.assertEqual(expected_participants, set(lang_O2_degradation_rxn.participants))

    def test_ozone_rxn(self):
        lang_Ozone_rxn = self.submodel.reactions.get_one(id='Ozone')
        self.assertEqual(lang_Ozone_rxn.reversible, True)
        self.assertEqual(len(lang_Ozone_rxn.rate_laws), 2)
        self.assertEqual(set(('1+2', '1+3')),
            set([rl.equation.expression for rl in lang_Ozone_rxn.rate_laws]))

    def test_rxn_participants(self):
        for lang_rxn in self.submodel.reactions:
            kb_rxn = self.kb.cell.reactions.get_one(id=lang_rxn.id)
            expected_participants = set()
            for part in kb_rxn.participants:
                lang_species = self.get_species(part.species.species_type.id, part.species.compartment.id)
                expected_participants.add(lang_species.species_coefficients.get_one(
                    coefficient=part.coefficient))
            self.assertEqual(expected_participants, set(lang_rxn.participants))

    def test_oxygen_ionization_rate_law(self):
        lang_Oxygen_ionization_rxn = self.submodel.reactions.get_one(id='Oxygen_ionization')
        lang_Oxygen_ionization_rate_law = lang_Oxygen_ionization_rxn.rate_laws[0]
        self.assertEqual(lang_Oxygen_ionization_rate_law.reaction, lang_Oxygen_ionization_rxn)
        self.assertEqual(lang_Oxygen_ionization_rate_law.direction, wc_lang.RateLawDirection.forward)
        self.assertEqual(lang_Oxygen_ionization_rate_law.k_cat, 1)
        self.assertEqual(lang_Oxygen_ionization_rate_law.k_m, 2)
        self.assertEqual(lang_Oxygen_ionization_rate_law.comments, 'custom')
        equation = lang_Oxygen_ionization_rate_law.equation
        self.assertEqual(equation.expression, '2 * O[c]')
        self.assertEqual(equation.parameters, [])
        lang_species = self.get_species('O', 'c')
        self.assertEqual(equation.modifiers, [lang_species])
