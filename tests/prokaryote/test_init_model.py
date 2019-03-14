from test.support import EnvironmentVarGuard
from wc_model_gen import prokaryote
from wc_utils.util.units import unit_registry
import unittest
import wc_kb
import wc_lang

class InitalizeModelTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        env = EnvironmentVarGuard()
        env.set('CONFIG__DOT__wc_kb__DOT__io__DOT__strict', '0')
        with env:
            cls.kb = wc_kb.io.Reader().run('tests/fixtures/test_broken_kb.xlsx',
                                           'tests/fixtures/test_broken_seq.fna',
                                            )[wc_kb.KnowledgeBase][0]

        cls.model = prokaryote.ProkaryoteModelGenerator(knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel]).run()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_compartments(self):
        cytosol = self.model.compartments.get_one(id='c')
        extracellular_space = self.model.compartments.get_one(id='e')
        self.assertIsInstance(cytosol, wc_lang.Compartment)
        self.assertIsInstance(extracellular_space, wc_lang.Compartment)

    def test_parameters(self):
        mean_doubling_time = self.model.parameters.get_one(id='mean_doubling_time')
        self.assertIsInstance(mean_doubling_time, wc_lang.Parameter)

    def test_metabolite_species(self):
        cytosol = self.model.compartments.get_one(id='c')
        for species in self.kb.cell.species_types.get(__type=wc_kb.core.MetaboliteSpeciesType):
            model_species_type = self.model.species_types.get_one(id=species.id)
            model_specie       = model_species_type.species.get_one(compartment=cytosol)

            self.assertIsInstance(model_species_type, wc_lang.SpeciesType)
            self.assertIsInstance(model_specie, wc_lang.Species, "Model does not contain species {}[c]".format(species.id))

    def test_rna_species(self):
        cytosol = self.model.compartments.get_one(id='c')
        for species in self.kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType):
            model_species_type = self.model.species_types.get_one(id=species.id)
            model_specie       = model_species_type.species.get_one(compartment=cytosol)

            self.assertIsInstance(model_species_type, wc_lang.SpeciesType)
            self.assertIsInstance(model_specie, wc_lang.Species)
            self.assertTrue(len(model_species_type.species)!=0)
            self.assertTrue(model_species_type.molecular_weight!=0)

    def test_protein_species(self):
        cytosol = self.model.compartments.get_one(id='c')
        for species in self.kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType):
            model_species_type = self.model.species_types.get_one(id=species.id)
            model_specie       = model_species_type.species.get_one(compartment=cytosol)

            self.assertIsInstance(model_species_type, wc_lang.SpeciesType)
            self.assertIsInstance(model_specie, wc_lang.Species)
            self.assertTrue(len(model_species_type.species)!=0)
            self.assertTrue(model_species_type.molecular_weight!=0)

    def test_concentrations(self):
        cytosol = self.model.compartments.get_one(id='c')
        for conc in self.kb.cell.concentrations:
            model_species_type = self.model.species_types.get_one(id=conc.species.species_type.id)
            model_specie = model_species_type.species.get_one(compartment=cytosol)

            self.assertIsInstance(model_species_type, wc_lang.SpeciesType)
            self.assertIsInstance(model_specie, wc_lang.Species)
            self.assertTrue(conc.units, unit_registry.parse_units('M'))
