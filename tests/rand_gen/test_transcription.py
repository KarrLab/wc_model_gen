""" Tests of transcription submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb_gen import random
from wc_model_gen.rand_gen import transcription, metabolism
import numpy
import scipy
import unittest
import wc_kb
import wc_lang


class TranscriptionSubmodelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = random.RandomKbGenerator(options={
            'component': {
                'PropertiesGenerator': {
                    'mean_volume': 1e-15,
                    'mean_doubling_time': 1000.,
                },
                'GenomeGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 100.,
                    'mean_gene_len': 10.,
                    'mean_copy_number': 10.,
                    'mean_half_life': 120.,
                },
                'MetabolitesGenerator': {
                },
            },
        }).run()
        cell = kb.cell
        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)

        model = wc_lang.Model()
        met = metabolism.MetabolismSubmodelGenerator(kb, model, options={})
        met.run()
        gen = transcription.TranscriptionSubmodelGenerator(
            kb, model, options={})
        gen.run()

        submodel = model.submodels.get_one(id='transcription')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')

        # check species types and species generated
        atp = model.species_types.get_one(id='atp')
        atp_cytosol = atp.species.get_one(compartment=cytosol)
        self.assertEqual(atp_cytosol.concentration.units,
                         wc_lang.ConcentrationUnit.M)

        concs = []
        for species_type in model.species_types:
            if species_type.id.startswith('rna_'):
                species = species_type.species.get_one(compartment=cytosol)
                concs.append(species.concentration.value)
                self.assertEqual(species.concentration.units,
                                 wc_lang.ConcentrationUnit.M)
        numpy.testing.assert_almost_equal(numpy.mean(
            concs), 10. / scipy.constants.Avogadro / 1e-15, decimal=2)

        # check reactions generated
        self.assertEqual(len(submodel.reactions), len(rnas))
        atp = model.species_types.get_one(
            id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(
            id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(
            id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(
            id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(
            id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(
            id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(
            id='h').species.get_one(compartment=cytosol)
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=atp).coefficient
            + submodel.reactions[0].participants.get_one(species=ctp).coefficient
            + submodel.reactions[0].participants.get_one(species=gtp).coefficient
            + submodel.reactions[0].participants.get_one(species=utp).coefficient,
            -rnas[0].get_len())
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=ppi).coefficient,
            rnas[0].get_len())
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h2o).coefficient,
            rnas[0].get_len() - 1)
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h).coefficient,
            -(rnas[0].get_len() - 1))

        rna_poly = model.observables.get_one(
            id='rna_poly_obs')
        rna_poly = rna_poly.species[0].species
        poly_avg_conc = 3000/scipy.constants.Avogadro / cytosol.initial_volume #http://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=2&id=106199
        # check rate laws
        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertEqual(rl.direction.name, 'forward')
            self.assertEqual(rl.equation.expression, '(((k_cat * rna_poly[c]) / (k_m + rna_poly[c])))')
            self.assertEqual(rl.equation.modifiers, [rna_poly])
            self.assertEqual(rl.equation.parameters, [])
            numpy.testing.assert_equal(rl.k_m, poly_avg_conc)
