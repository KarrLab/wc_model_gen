""" Tests of metabolism submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb_gen
import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_lang
import wc_kb


class MetabolismSubmodelGeneratorTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/min_kb.xlsx',
                                       'tests/fixtures/min_kb_seq.fna',
                                        strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(
                        knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel,
                                              prokaryote.MetabolismSubmodelGenerator],
                        options = {'component': {
                             'TranscriptionSubmodelGenerator': {
                               'rate_dynamics': 'phenomenological'}}}).run()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_metabolite_species(self):

        cytosol = self.model.compartments.get_one(id='c')
        for species in self.kb.cell.species_types.get(__type=wc_kb.core.MetaboliteSpeciesType):
            model_species_type = self.model.species_types.get_one(id=species.id)
            model_specie       = model_species_type.species.get_one(compartment=cytosol)

            self.assertIsInstance(model_species_type, wc_lang.SpeciesType)
            self.assertIsInstance(model_specie, wc_lang.Species)
