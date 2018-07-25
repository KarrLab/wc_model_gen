""" Tests of RNA degradation submodel generation

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-24
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb_gen import random
from wc_model_gen.rand_gen import protein_degradation, metabolism
import numpy
import scipy
import unittest
import wc_kb
import wc_lang


class ProteinDegradationSubmodelGeneratorTestCase(unittest.TestCase):
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

        prots = cell.species_types.get(__type=wc_kb.ProteinSpeciesType)

        model = wc_lang.Model()
        met = metabolism.MetabolismSubmodelGenerator(kb, model, options={})
        met.run()
        gen = protein_degradation.ProteinDegradationSubmodelGenerator(
            kb, model, options={})
        gen.run()

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
        self.assertEqual(submodel.reactions[0].participants.get_one(species=atp).coefficient, -1)
        self.assertEqual(submodel.reactions[0].participants.get_one(species=adp).coefficient, 1)
        self.assertEqual(submodel.reactions[0].participants.get_one(species=pi).coefficient, 1)
        self.assertEqual(
            + submodel.reactions[0].participants.get_one(species=h2o).coefficient,
            -(prots[0].get_len()-1))

        aa_species = model.species_types.get_one(
                    id='cys').species.get_one(compartment=cytosol)
        self.assertEqual(submodel.reactions[0].participants.get_one(species = aa_species).coefficient, prots[0].get_seq().count('C'))
        

        # check rate laws
        deg_protease = model.observables.get_one(id='deg_protease_obs')
        deg_protease = deg_protease.species[0].species.species_type
        deg_avg_conc = 5000/scipy.constants.Avogadro / cytosol.initial_volume
        for prot, rxn in zip(prots, submodel.reactions):
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertEqual(rl.direction.name, 'forward')
            self.assertEqual(rl.equation.expression,
                '{0}[c] * (((k_cat * {1}[c]) / (k_m + {1}[c])) + {2})'.format(prot.id, deg_protease.id, '0.1'))
            if prot.id != deg_protease.id:
                self.assertEqual(rl.equation.modifiers, [deg_protease.species.get_one(compartment=cytosol), rxn.participants[0].species])
            else:
                self.assertEqual(rl.equation.modifiers, [rxn.participants[0].species])
            self.assertEqual(rl.equation.parameters, [])
            self.assertEqual(rl.k_m,deg_avg_conc)
            self.assertEqual(rl.k_cat, 2 * numpy.log(2) / prot.half_life)

