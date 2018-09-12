""" Testing Translation Submodel Generator

:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-07-23
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb_gen
import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_lang
import wc_kb


class TranslationSubmodelGeneratorTestCase(unittest.TestCase):

    def test(self):
        kb = wc_kb_gen.random.RandomKbGenerator(options={
            'component': {
                'PropertiesGenerator': {
                    'mean_volume': 1e-15,
                    'mean_cell_cycle_length': 1000.,
                },
                'GenomeGenerator': {
                    'num_chromosomes': 1,
                    'mean_num_genes': 200,
                    'mean_gene_len': 10,
                    'mean_copy_number': 10,
                    'mean_half_life': 120,
                },
                'MetabolitesGenerator': {
                },
            },
        }).run()
        cell = kb.cell

        model = prokaryote.ProkaryoteModelGenerator(
                     knowledge_base=kb,
                     component_generators=[prokaryote.InitalizeModel,
                                           prokaryote.TranslationSubmodelGenerator]).run()

        submodel = model.submodels.get_one(id='translation')

        # check compartments generated
        cytosol = model.compartments.get_one(id='c')
        self.assertEqual(cytosol.name, 'Cytosol')

        for species_type in model.species_types:
            if species_type.id.startswith('prot_'):
                species = species_type.species.get_one(compartment=cytosol)
                self.assertEqual(species.concentration.units, wc_lang.ConcentrationUnit.M)

        # check reactions generated
        prots = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        self.assertEqual(len(submodel.reactions), len(prots))

        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        gdp = model.species_types.get_one(id='gdp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)
        ribosome = model.observables.get_one(id='complex_70S_obs').expression.species[0]
        initiation_factors = model.observables.get_one(id='translation_init_factors_obs').expression.species[0]
        elongation_factors = model.observables.get_one(id='translation_elongation_factors_obs').expression.species[0]
        release_factors = model.observables.get_one(id='translation_release_factors_obs').expression.species[0]

        for reaction in submodel.reactions:
            prot_model = model.species_types.get_one(id=reaction.id[12:])
            prot_kb = kb.cell.species_types.get_one(id=reaction.id[12:])
            length = len(prot_kb.get_seq())

            self.assertEqual(reaction.participants.get_one(species=prot_model.species[0]).coefficient, 1)
            self.assertEqual(reaction.participants.get_one(species=gtp).coefficient, -(length+2))
            self.assertEqual(reaction.participants.get_one(species=gdp).coefficient, (length+2))
            self.assertEqual(reaction.participants.get_one(species=pi).coefficient, 2*length)

            # TODO: following species are on both sides of reaction, thus gettign error
            #self.assertEqual(reaction.participants.get_one(species=initiation_factors).coefficient, 1)
            #self.assertEqual(reaction.participants.get_one(species=elongation_factors).coefficient, length)
            #self.assertEqual(reaction.participants.get_one(species=release_factors).coefficient, 1)

    @unittest.skip
    def test(self):
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
