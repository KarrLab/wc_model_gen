""" Testing Translation Submodel Generator

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-23
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_kb_gen import random
from wc_model_gen.prokaryote import translation, metabolism
import numpy
import scipy
import unittest
import wc_kb
import wc_lang


class TranslationSubmodelGeneratorTestCase(unittest.TestCase):
    def test(self):
        kb = random.RandomKbGenerator(options={
            'component': {
                'PropertiesGenerator': {
                    'mean_volume': 1e-15,
                    'mean_doubling_time': 1000.,
                },
                'GenomeGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 200.,
                    'mean_gene_len': 10.,
                    'mean_copy_number': 10.,
                    'mean_half_life': 120.,
                },
                'MetabolitesGenerator': {
                },
            },
        }).run()
        cell = kb.cell

        model = wc_lang.Model()
        met = metabolism.MetabolismSubmodelGenerator(kb, model, options={})
        met.run()
        gen = translation.TranslationSubmodelGenerator(kb, model, options={})
        gen.run()

        submodel = model.submodels.get_one(id='translation')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'cytosol')

        for species_type in model.species_types:
            if species_type.id.startswith('prot_'):
                species = species_type.species.get_one(compartment=cytosol)
                self.assertEqual(species.concentration.units,
                                 wc_lang.ConcentrationUnit.M)

        # check reactions generated
        prots = cell.species_types.get(__type=wc_kb.ProteinSpeciesType)

        self.assertEqual(len(submodel.reactions), 3 * len(prots))
        gtp = model.species_types.get_one(
            id='gtp').species.get_one(compartment=cytosol)
        gdp = model.species_types.get_one(
            id='gdp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(
            id='pi').species.get_one(compartment=cytosol)
        complex_70S = model.species_types.get_one(
            id='complex_70S_A').species.get_one(compartment=cytosol)
        complex_70S_I = model.species_types.get_one(
            id='complex_70S_IA').species.get_one(compartment=cytosol)

        for reaction in submodel.reactions:
            prot = cell.species_types.get_one(id=reaction.name)
            length = len(prot.get_seq())

            if reaction.id.startswith('translation_init_'):  # translation init
                self.assertEqual(reaction.participants.get_one(
                    species=gtp).coefficient, -1)
                self.assertEqual(reaction.participants.get_one(
                    species=gdp).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(
                    species=pi).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(
                    species=complex_70S).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(
                    species=complex_70S_I).coefficient, -1)

            elif reaction.id.startswith('translation_elon_'):  # translation elon
                self.assertEqual(reaction.participants.get_one(
                    species=gtp).coefficient, -2 * length)
                self.assertEqual(reaction.participants.get_one(
                    species=gdp).coefficient, 2 * length)
                self.assertEqual(reaction.participants.get_one(
                    species=pi).coefficient, 2 * length)
                for participant in reaction.participants:
                    if participant.species.species_type.id == prot.id+'_att':
                        self.assertEqual(participant.coefficient, 1)
                        break

            else:  # translation term
                self.assertEqual(reaction.participants.get_one(
                    species=gtp).coefficient, -1)
                self.assertEqual(reaction.participants.get_one(
                    species=gdp).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(
                    species=pi).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(
                    species=complex_70S_I).coefficient, 1)
                for participant in reaction.participants:
                    if participant.species.species_type.id == prot.id+'_att':
                        self.assertEqual(participant.coefficient, -1)

   """     # check rate laws
        for rxn in submodel.reactions:
            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            if rxn.id.startswith('translation_init_'):
                exp = 'k_cat * (IF_obs' + \
                      '/ (k_m +IF_obs))'
                self.assertEqual(rl.equation.observables, [
                                 model.observables.get_one(id='IF_obs')])
                prot_id = rxn.id[rxn.id.find('init_')+5:]
                prot = cell.species_types.get_one(id=prot_id)
                self.assertEqual(
                    rl.k_cat, 2 * (numpy.log(2) / prot.half_life + numpy.log(2) / 1000))
            elif rxn.id.startswith('translation_elon_'):
                exp = 'k_cat * (EF_obs' + \
                      '/ (k_m +EF_obs))'
                self.assertEqual(rl.equation.observables, [
                                 model.observables.get_one(id='EF_obs')])
                prot_id = rxn.id[rxn.id.find('elon_')+5:]
                prot = cell.species_types.get_one(id=prot_id)
                self.assertEqual(
                    rl.k_cat, 2 * (numpy.log(2) / prot.half_life + numpy.log(2) / 1000))
            else:
                exp = 'k_cat * (RF_obs' + \
                      '/ (k_m +RF_obs))'
                self.assertEqual(rl.equation.observables, [
                                 model.observables.get_one(id='RF_obs')])
                prot_id = rxn.id[rxn.id.find('term_')+5:]
                prot = cell.species_types.get_one(id=prot_id)
                self.assertEqual(
                    rl.k_cat, 2 * (numpy.log(2) / prot.half_life + numpy.log(2) / 1000))
            self.assertEqual(rl.direction.name, 'forward')
            self.assertEqual(rl.equation.expression, exp)
            self.assertEqual(rl.equation.modifiers, [])
            self.assertEqual(rl.k_m, 0.05)
"""
