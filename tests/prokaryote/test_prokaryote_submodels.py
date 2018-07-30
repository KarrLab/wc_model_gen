import os
import shutil
import unittest
import tempfile
import wc_kb
import wc_sim
import wc_test
import wc_lang
import wc_model_gen
import wc_model_gen.prokaryote as prokaryote

class ProkaryoteTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx','tests/fixtures/seq.fna', strict = False)
        cls.model = wc_model_gen.ModelGenerator(
                            knowledge_base=cls.kb,
                            component_generators=[prokaryote.CompartmentsGenerator,
                                                  prokaryote.ParametersGenerator,
                                                  prokaryote.MetaboliteSpeciesGenerator]).run()

        cls.model.id = 'prokaryote_submodels_tests'
        cls._results_dir = tempfile.mkdtemp()
        cls.model.version = '0.0.1'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._results_dir)

    def test_compartments_generator(self):
        self.assertIsInstance(self.model.compartments.get_one(id='c'), wc_lang.core.Compartment)
        self.assertIsInstance(self.model.compartments.get_one(id='e'), wc_lang.core.Compartment)

    def test_parameters_generator(self):
        self.assertIsInstance(self.model.parameters.get_one(id='cellCycleLength'), wc_lang.core.Parameter)
        self.assertIsInstance(self.model.parameters.get_one(id='fractionDryWeight'), wc_lang.core.Parameter)

    def test_metabolites_generator(self):
        comp = self.model.compartments.get_one(id='c')
        specie = self.model.species_types.get_one(id='ATP').species.get_one(compartment=comp)
        self.assertIsInstance(specie, wc_lang.core.Species)

        comp = self.model.compartments.get_one(id='c')
        specie = self.model.species_types.get_one(id='Ala').species.get_one(compartment=comp)
        self.assertIsInstance(specie, wc_lang.core.Species)

        comp = self.model.compartments.get_one(id='c')
        specie = self.model.species_types.get_one(id='H2O').species.get_one(compartment=comp)
        self.assertIsInstance(specie, wc_lang.core.Species)

        comp = self.model.compartments.get_one(id='e')
        specie = self.model.species_types.get_one(id='H2O').species.get_one(compartment=comp)
        self.assertIsInstance(specie, wc_lang.core.Species)

    def test_construction(self):
        self.assertEqual(self.model.get_reactions(), [])
        self.assertEqual(self.model.get_rate_laws(), [])
        self.assertEqual(self.model.submodels, [])

class TranscriptionSubmodelTestCase(ProkaryoteTestCase):

    @classmethod
    def setUpClass(cls):
        super(TranscriptionSubmodelTestCase, cls).setUpClass()
        prokaryote.TranscriptionSubmodelGenerator(cls.kb, cls.model).run()

    def test_construction(self):
        self.assertEqual(
            len(self.model.get_species_types(type=wc_lang.SpeciesTypeType.rna)),
            len(self.kb.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType)))

        self.assertEqual(
            len(self.model.get_reactions()),
            len(self.kb.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType)))

        self.assertEqual(
            len(self.model.get_rate_laws()),
            len(self.kb.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType)))

class TranslationSubmodelTestCase(ProkaryoteTestCase):

    @classmethod
    def setUpClass(cls):
        super(TranslationSubmodelTestCase, cls).setUpClass()
        prokaryote.TranslationSubmodelGenerator(cls.kb, cls.model).run()

    def test_construction(self):
        self.assertEqual(
            len(self.model.get_species_types(type=wc_lang.SpeciesTypeType.protein)),
            len(self.kb.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType)))

        self.assertEqual(
            len(self.model.get_reactions()),
            3*len(self.kb.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType)))

        self.assertEqual(
            len(self.model.get_rate_laws()),
            3*len(self.kb.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType)))

class DegradationSubmodelTestCase(ProkaryoteTestCase):

    @classmethod
    def setUpClass(cls):
        super(DegradationSubmodelTestCase, cls).setUpClass()
        prokaryote.DegradationSubmodelGenerator(cls.kb, cls.model).run()

    def test_construction(self):
        # These tests will be made more meaningfull after finishign the ability to reduce KB and hence model size
        self.assertEqual(len(self.model.get_species_types()), 37)
        self.assertEqual(len(self.model.get_reactions()), 0)
        self.assertEqual(len(self.model.get_rate_laws()), 0)
