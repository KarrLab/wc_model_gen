""" Tests of RNA degradation submodel generation

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import rand_wc_model_gen
from rand_wc_model_gen import kb_gen
from wc_model_gen.rand_gen import rna_degradation, metabolism
import numpy
import scipy
import unittest
import wc_kb
import wc_lang


class RnaDegradationSubmodelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = kb_gen.KbGenerator(options={
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
        gen = rna_degradation.RnaDegradationSubmodelGenerator(
            kb, model, options={})
        gen.run()

        submodel = model.submodels.get_one(id='rna_degradation')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')

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
        numpy.testing.assert_almost_equal(numpy.mean(
            concs), 10. / scipy.constants.Avogadro / 1e-15, decimal=2)

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

        # check rate laws
        '''for rna, rxn in zip(rnas, submodel.reactions):
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertEqual(rl.direction.name, 'forward')
            self.assertEqual(rl.equation.expression,
                             'k_cat * {0}[c] / (k_m + {0}[c])'.format(rna.id))
            self.assertEqual(rl.equation.modifiers, [
                             rxn.participants[0].species])
            self.assertEqual(rl.equation.parameters, [])
            self.assertEqual(rl.k_m, rna.concentration)

        k_cats = [rxn.rate_laws[0].k_cat for rxn in submodel.reactions]
        numpy.testing.assert_almost_equal(numpy.mean(
            k_cats), 2. * numpy.log(2.) * (1. / 120.), decimal=2)'''
