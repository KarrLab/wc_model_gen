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

    def test(self):
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
                                           prokaryote.MetabolismSubmodelGenerator]).run()

        cell = kb.cell
        prots = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        submodel = model.submodels.get_one(id='metabolism')

        # check parameters generated
        self.assertEqual(model.parameters.get_one(id='fraction_dry_weight').value, 0.3)

        # check species types and species generated
        cytosol = model.compartments.get(id='c')[0]

        for species in kb.cell.species_types.get(__type=wc_kb.core.MetaboliteSpeciesType):
            model_species = model.species_types.get_one(id=species.id)
            model_species_cytosol = model_species.species.get_one(
                compartment=cytosol)
            self.assertIsInstance(model_species, wc_lang.SpeciesType)
            self.assertIsInstance(model_species_cytosol, wc_lang.Species)
