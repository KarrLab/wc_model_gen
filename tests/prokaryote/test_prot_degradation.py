""" Tests of RNA degradation submodel generation

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-24
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb_gen
import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_lang
import wc_kb

class ProteinDegradationSubmodelGeneratorTestCase(unittest.TestCase):

    def test(self):
        """
        kb = wc_kb.io.Reader().run('tests/fixtures/core.xlsx',
                                            'tests/fixtures/seq.fna', strict=False)

        model = prokaryote.ProkaryoteModelGenerator(
                    knowledge_base = kb,
                    component_generators = [prokaryote.InitalizeModel,
                                            prokaryote.ProteinDegradationSubmodelGenerator]).run()
        """

        kb = wc_kb_gen.random.RandomKbGenerator(options={
             'component': {
                 'GenomeGenerator': {
                     'mean_num_genes': 20,
                     'mean_gene_len': 50,
                     'num_ncRNA': 0,
                     'translation_table': 4,
                     'mean_copy_number': 100,
                     'mean_half_life': 100
                 },
                 'PropertiesGenerator': {
                     'mean_cell_cycle_length': 100,
                 },
             },
         }).run()

        model = prokaryote.ProkaryoteModelGenerator(
                     knowledge_base=kb,
                     component_generators=[prokaryote.InitalizeModel,
                                           prokaryote.ProteinDegradationSubmodelGenerator]).run()

        cell = kb.cell
        prots = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        submodel = model.submodels.get_one(id='protein_degradation')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'Cytosol')

        # check species types and species generated        
        for species in kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType):
            model_species = model.species_types.get_one(id=species.id)
            model_species_cytosol = model_species.species.get_one(
                compartment=cytosol)
            self.assertIsInstance(model_species, wc_lang.SpeciesType)
            self.assertIsInstance(model_species_cytosol, wc_lang.Species)

        # check reactions generated
        self.assertEqual(len(submodel.reactions), len(prots))
        atp = model.species_types.get_one(
            id='atp').species.get_one(compartment=cytosol)
        adp = model.species_types.get_one(
            id='adp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(
            id='pi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(
            id='h2o').species.get_one(compartment=cytosol)
        self.assertEqual(submodel.reactions[0].participants.get_one(
            species=atp).coefficient, -1)
        self.assertEqual(submodel.reactions[0].participants.get_one(
            species=adp).coefficient, 1)
        self.assertEqual(submodel.reactions[0].participants.get_one(
            species=pi).coefficient, 1)
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h2o).coefficient,
            -(prots[0].get_len()-1))

        aa_species = model.species_types.get_one(
            id='cys').species.get_one(compartment=cytosol)
        self.assertEqual(submodel.reactions[0].participants.get_one(
            species=aa_species).coefficient, prots[0].get_seq().count('C'))

        for prot, rxn in zip(prots, submodel.reactions):
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertEqual(rl.direction.name, 'forward')
            self.assertIsInstance(rl.k_m, float)
            self.assertIsInstance(rl.k_cat, float)
