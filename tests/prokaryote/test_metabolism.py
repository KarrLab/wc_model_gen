""" Tests of metabolism submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from test.support import EnvironmentVarGuard
import wc_kb_gen
from wc_model_gen import prokaryote
import unittest
import wc_lang
import wc_kb

class MetabolismSubmodelGeneratorTestCase(unittest.TestCase):

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
                                              prokaryote.TranscriptionSubmodelGenerator,
                                              prokaryote.MetabolismSubmodelGenerator]).run()

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

    def test_raises(self):

        break_model = self.model
        rxn = break_model.submodels.get_one(id='metabolism').reactions[0]
        rxn.id = 'wrong_id_name'
        submod_gen = prokaryote.MetabolismSubmodelGenerator(
                        knowledge_base = self.kb,
                        model = self.model)
        self.assertRaisesRegexp(Exception, 'invalid reaction id, no associated rate law is defined', submod_gen.gen_rate_laws)

        break_model = self.model
        rxn = break_model.submodels.get_one(id='metabolism').reactions[0]
        rxn.id = 'transfer_wrong_id_name'
        submod_gen = prokaryote.MetabolismSubmodelGenerator(
                        knowledge_base = self.kb,
                        model = self.model)
        self.assertRaisesRegexp(Exception, 'invalid reaction id, no associated rate law is defined', submod_gen.gen_rate_laws)

        break_model = self.model
        break_model.species_types.get_one(id='atp').id='pta'
        submod_gen = prokaryote.MetabolismSubmodelGenerator(
                        knowledge_base = self.kb,
                        model = self.model)
        self.assertRaisesRegexp(ValueError, 'not found in model.species', submod_gen.gen_reactions)
