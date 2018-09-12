""" Tests of RNA degradation submodel generation

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-24
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb_gen
from wc_model_gen.prokaryote import protein_degradation, metabolism
import wc_model_gen.prokaryote as prokaryote
import numpy
import scipy
import unittest
import wc_kb
import wc_lang


class ProteinDegradationSubmodelGeneratorTestCase(unittest.TestCase):
    def test(self):

        rand_kb = wc_kb_gen.random.RandomKbGenerator(options={
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
                     'mean_doubling_time': 100,
                 },
             },
         }).run()

        cell = rand_kb.cell
        prots = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)

        model = prokaryote.ProkaryoteModelGenerator(
                     knowledge_base=rand_kb,
                     component_generators=[prokaryote.CompartmentsGenerator,
                                           prokaryote.ParametersGenerator,
                                           prokaryote.MetabolismSubmodelGenerator,
                                           prokaryote.TranscriptionSubmodelGenerator,
                                           prokaryote.TranslationSubmodelGenerator,
                                           prokaryote.ProteinDegradationSubmodelGenerator]).run()

        submodel = model.submodels.get_one(id='protein_degradation')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')

        # check species types and species generated
        atp = model.species_types.get_one(id='atp')
        atp_cytosol = atp.species.get_one(compartment=cytosol)
        self.assertEqual(atp_cytosol.concentration.units,
                         wc_lang.ConcentrationUnit.M)

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

        # check rate laws
        deg_protease = model.observables.get_one(id='deg_protease_obs')
        deg_avg_conc = 5000/scipy.constants.Avogadro / cytosol.initial_volume

        for prot, rxn in zip(prots, submodel.reactions):
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertEqual(rl.direction.name, 'forward')
            self.assertIsInstance(rl.k_m, float)
            self.assertIsInstance(rl.k_cat, float)
