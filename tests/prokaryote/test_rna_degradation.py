""" Tests of RNA degradation submodel generation

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

class RnaDegradationSubmodelGeneratorTestCase(unittest.TestCase):
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
                                           prokaryote.RnaDegradationSubmodelGenerator]).run()


        cell = kb.cell
        rnas = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)

        submodel = model.submodels.get_one(id='rna_degradation')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'Cytosol')

        # check species types and species generated
        amp = model.species_types.get_one(id='amp')
        amp_cytosol = amp.species.get_one(compartment=cytosol)
        self.assertEqual(amp_cytosol.concentration.units,
                         wc_lang.ConcentrationUnit.M)

        concs = []
        for species_type in model.species_types:
            if species_type.id.startswith('rna_'):
                species = species_type.species.get_one(compartment=cytosol)
                concs.append(species.concentration.value)
                self.assertEqual(species.concentration.units,
                                 wc_lang.ConcentrationUnit.M)

        # check reactions generated
        self.assertEqual(len(submodel.reactions), len(rnas))
        amp = model.species_types.get_one(
            id='amp').species.get_one(compartment=cytosol)
        cmp = model.species_types.get_one(
            id='cmp').species.get_one(compartment=cytosol)
        gmp = model.species_types.get_one(
            id='gmp').species.get_one(compartment=cytosol)
        ump = model.species_types.get_one(
            id='ump').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(
            id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(
            id='h').species.get_one(compartment=cytosol)
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=amp).coefficient
            + submodel.reactions[0].participants.get_one(species=cmp).coefficient
            + submodel.reactions[0].participants.get_one(species=gmp).coefficient
            + submodel.reactions[0].participants.get_one(species=ump).coefficient,
            rnas[0].get_len())
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h2o).coefficient,
            -(rnas[0].get_len() - 1))
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h).coefficient,
            rnas[0].get_len() - 1)

    """ # check rate laws
     deg_rnase = model.observables.get_one(id='deg_rnase_obs')
     deg_avg_conc = 5000/scipy.constants.Avogadro / cytosol.initial_volume
     for rna, rxn in zip(rnas, submodel.reactions):
         self.assertEqual(len(rxn.rate_laws), 1)
         rl = rxn.rate_laws[0]
         self.assertEqual(rl.direction.name, 'forward')
         self.assertEqual(rl.k_m, deg_avg_conc)
         self.assertEqual(rl.k_cat, 2 * numpy.log(2) / rna.half_life)"""
