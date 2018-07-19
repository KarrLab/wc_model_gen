


import rand_wc_model_gen
from rand_wc_model_gen import kb_gen
from wc_model_gen.rand_gen import translation, metabolism
import numpy
import scipy
import unittest
import wc_kb
import wc_lang


class TranslationSubmodelGeneratorTestCase(unittest.TestCase):   
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
                self.assertEqual(species.concentration.units, wc_lang.ConcentrationUnit.M)

        # check reactions generated
        prots = cell.species_types.get(__type=wc_kb.ProteinSpeciesType)
        
        self.assertEqual(len(submodel.reactions), 3 * len(prots))
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        gdp = model.species_types.get_one(id='gdp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)
        complex_70S = model.species_types.get_one(id='complex_70S_A').species.get_one(compartment=cytosol)
        complex_70S_I = model.species_types.get_one(id='complex_70S_IA').species.get_one(compartment=cytosol)


        for reaction in submodel.reactions:
            prot = cell.species_types.get_one(id=reaction.name)
            length = len(prot.get_seq())

            if reaction.id.startswith('translation_init_'):     #translation init
                self.assertEqual(reaction.participants.get_one(species = gtp).coefficient, -1)
                self.assertEqual(reaction.participants.get_one(species = gdp).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(species = pi).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(species = complex_70S).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(species = complex_70S_I).coefficient, -1)

            elif reaction.id.startswith('translation_elon_'): #translation elon
                self.assertEqual(reaction.participants.get_one(species = gtp).coefficient, -2 * length)
                self.assertEqual(reaction.participants.get_one(species = gdp).coefficient, 2 * length)
                self.assertEqual(reaction.participants.get_one(species = pi).coefficient, 2 * length)
                for participant in reaction.participants:
                    if participant.species.species_type.id == prot.id+'_att':
                        self.assertEqual(participant.coefficient, 1)
                        break

            else: #translation term
                self.assertEqual(reaction.participants.get_one(species = gtp).coefficient, -1)
                self.assertEqual(reaction.participants.get_one(species = gdp).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(species = pi).coefficient, 1)
                self.assertEqual(reaction.participants.get_one(species = complex_70S_I).coefficient, 1)
                for participant in reaction.participants:
                    if participant.species.species_type.id == prot.id+'_att':
                        self.assertEqual(participant.coefficient, -1)

                
        # check rate laws
        '''for rxn in submodel.reactions:
            exp = 'k_cat'
            for participant in rxn.participants:
                if participant.coefficient < 0:
                    exp = exp + ' * (' + participant.species.id() + '/ (k_m + ' + participant.species.id() + '))'

            self.assertEqual(len(rxn.rate_laws), 1)
            rl = rxn.rate_laws[0]
            self.assertEqual(rl.direction.name, 'forward')
            self.assertEqual(rl.equation.expression, exp)
            self.assertEqual(rl.equation.parameters, [])            
            self.assertEqual(rl.k_m, 1)
            self.assertEqual(rl.k_cat, 1)'''
