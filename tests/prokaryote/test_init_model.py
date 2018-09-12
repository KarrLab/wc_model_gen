import shutil
import os
import unittest
import tempfile
import wc_kb
import wc_lang
import wc_model_gen
import wc_model_gen.prokaryote as prokaryote


class InitalizeModelTestCase(unittest.TestCase):

    @unittest.skip(' coderefractoring')
    def test_initalize_model(self):

        """
        def setUp(self):
            self.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx',
                                            'tests/fixtures/seq.fna', strict=False)
            self.model = wc_model_gen.ModelGenerator(
                knowledge_base=self.kb,
                component_generators=[prokaryote.InitalizeModel]).run()

            self.model.id = 'prokaryote_compartment_tests'
            self._results_dir = tempfile.mkdtemp()

        def tearDown(self):
            shutil.rmtree(self._results_dir)

        def test_compartments(self):
            cytosol = self.model.compartments.get_one(id='c')
            self.assertEqual(cytosol.name, 'cytosol')
            ec = self.model.compartments.get_one(id='e')
            self.assertEqual(ec.name, 'extracellular space')


            class ParametersGeneratorTestCase(unittest.TestCase):

                def setUp(self):
                    self.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx',
                                                    'tests/fixtures/seq.fna', strict=False)
                    self.model = wc_model_gen.ModelGenerator(
                        knowledge_base=self.kb,
                        component_generators=[prokaryote.ParametersGenerator]).run()

                    self.model.id = 'prokaryote_compartment_tests'
                    self._results_dir = tempfile.mkdtemp()

                def tearDown(self):
                    shutil.rmtree(self._results_dir)

                def test_parameters(self):
                    self.assertIsInstance(self.model.parameters.get_one(
                        id='cellCycleLength'), wc_lang.core.Parameter)
                    self.assertIsInstance(self.model.parameters.get_one(
                        id='fractionDryWeight'), wc_lang.core.Parameter)





            class InitalizeModelTestCase(unittest.TestCase):

                def setUpClass(self):

                    self.kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx',
                                                    'tests/fixtures/seq.fna', strict=False)

                    self.model = prokaryote.ProkaryoteModelGenerator(
                                 knowledge_base=self.kb,
                                 component_generators=[prokaryote.InitalizeModel]).run()

                def testSpecies(self):

                            cytosol = self.model.compartments.get_one(id='c')
                            self.assertEqual(cytosol.name, 'cytosol')
                            ec = self.model.compartments.get_one(id='e')
                            self.assertEqual(ec.name, 'extracellular space')

                    self.assertIsInstance(self.model.parameters.get_one(
                        id='cellCycleLength'), wc_lang.core.Parameter)
                    self.assertIsInstance(self.model.parameters.get_one(
                        id='fractionDryWeight'), wc_lang.core.Parameter)


                    # check parameters generated
                    self.assertEqual(self.model.parameters.get_one(
                        id='fraction_dry_weight').value, 0.7)

                    # check species types and species generated
                    cytosol = self.model.compartments.get(id='c')[0]

                    for species in self.kb.cell.species_types.get(__type=wc_kb.core.MetaboliteSpeciesType):
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

                def test_rxn_participants(self):
                    for lang_rxn in self.submodel.reactions:
                        kb_rxn = self.kb.cell.reactions.get_one(id=lang_rxn.id)
                        expected_participants = set()
                        for part in kb_rxn.participants:
                            lang_species = self.get_species(part.species.species_type.id, part.species.compartment.id)
                            expected_participants.add(lang_species.species_coefficients.get_one(
                                coefficient=part.coefficient))
                        self.assertEqual(expected_participants, set(lang_rxn.participants))
            """
