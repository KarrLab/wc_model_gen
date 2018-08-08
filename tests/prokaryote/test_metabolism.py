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

    def testModelGen(self):
        submodel = self.model.submodels.get_one(id='metabolism')
        self.assertIsInstance(submodel, wc_lang.Submodel)

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
