import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_kb
import wc_lang

class InitalizeModelTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/test_broken.xlsx',
                                       'tests/fixtures/test_broken_seq.fna',
                                        strict=False)

        cls.model = prokaryote.ProkaryoteModelGenerator(knowledge_base = cls.kb,
                        component_generators=[prokaryote.InitalizeModel]).run()

    @classmethod
    def tearDownClass(cls):
        pass

    def test_compartments(self):
        cytosol = self.model.compartments.get_one(id='c')
        extracellular_space = self.model.compartments.get_one(id='e')
        self.assertIsInstance(cytosol, wc_lang.core.Compartment)
        self.assertIsInstance(extracellular_space, wc_lang.core.Compartment)

    def test_parameters(self):
        cell_cycle_length = self.model.parameters.get_one(id='cell_cycle_length')
        fraction_dry_weight = self.model.parameters.get_one(id='fraction_dry_weight')
        self.assertIsInstance(cell_cycle_length, wc_lang.core.Parameter)
        self.assertIsInstance(fraction_dry_weight, wc_lang.core.Parameter)

    def test_metabolite_species(self):
        cytosol = self.model.compartments.get_one(id='c')
        for species in self.kb.cell.species_types.get(__type=wc_kb.core.MetaboliteSpeciesType):
            model_species_type = self.model.species_types.get_one(id=species.id)
            model_specie       = model_species_type.species.get_one(compartment=cytosol)

            self.assertIsInstance(model_species_type, wc_lang.SpeciesType)
            self.assertIsInstance(model_specie, wc_lang.Species)

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
            self.assertTrue(conc.units, wc_lang.ConcentrationUnit.M)        
