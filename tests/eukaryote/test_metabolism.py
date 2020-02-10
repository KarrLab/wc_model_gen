""" Tests of metabolism submodel generation
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2020-02-03
:Copyright: 2019-2020, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import metabolism
from wc_onto import onto as wc_ontology
import collections
import conv_opt
import math
import unittest
import wc_lang
import wc_kb


class MetabolismSubmodelGeneratorTestCase(unittest.TestCase):

    def setUp(self):
        kb = self.kb = wc_kb.KnowledgeBase()
        model = self.model = wc_lang.Model()

        model.parameters.create(id='mean_doubling_time', value=math.log(2)/12)

        c = model.compartments.create(id='c')

        # Create metabolites
        for m in [1,2,3]:
            metabolite_st = model.species_types.create(id='m{}'.format(m))
            metabolite_species = model.species.create(species_type=metabolite_st, compartment=c)
            metabolite_species.id = metabolite_species.gen_id()

        # Create enzymes
        for i, conc in {'enzyme1':0.01, 'enzyme2': 0.01, 'enzyme3': 0.}.items():
            enzyme_st = model.species_types.create(id=i)
            enzyme_species = model.species.create(species_type=enzyme_st, compartment=c)
            enzyme_species.id = enzyme_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=enzyme_species, mean=conc)
            conc_model.id = conc_model.gen_id()

        # Create reactions in metabolism submodel
        self.gen = metabolism.MetabolismSubmodelGenerator(kb, model)        

        submodel = model.submodels.get_one(id='metabolism')
        
        ex1 = submodel.reactions.create(id='ex_m1', reversible=False, model=model)           
        ex1.participants.append(model.species.get_one(id='m1[c]').species_coefficients.get_or_create(coefficient=1))

        ex2 = submodel.reactions.create(id='ex_m2', reversible=False, model=model)           
        ex2.participants.append(model.species.get_one(id='m2[c]').species_coefficients.get_or_create(coefficient=1))

        ex3 = submodel.reactions.create(id='ex_m3', reversible=False, model=model)           
        ex3.participants.append(model.species.get_one(id='m3[c]').species_coefficients.get_or_create(coefficient=1))      
        
        r1 = submodel.reactions.create(id='r1', reversible=True, model=model)           
        r1.participants.append(model.species.get_one(id='m1[c]').species_coefficients.get_or_create(coefficient=-1))
        r1.participants.append(model.species.get_one(id='m2[c]').species_coefficients.get_or_create(coefficient=-1))
        r1.participants.append(model.species.get_one(id='m3[c]').species_coefficients.get_or_create(coefficient=1))
        r1_rate_law_expression1, error = wc_lang.RateLawExpression.deserialize('k_cat_r1_forward_enzyme3 * enzyme3[c]', {
            wc_lang.Parameter: {'k_cat_r1_forward_enzyme3': model.parameters.create(id='k_cat_r1_forward_enzyme3')},
            wc_lang.Species: {'enzyme3[c]': model.species.get_one(id='enzyme3[c]')},
            })
        assert error is None, str(error)
        r1_model_rate_law1 = model.rate_laws.create(
            expression=r1_rate_law_expression1,
            reaction=r1,
            direction=wc_lang.RateLawDirection['forward'])
        r1_model_rate_law1.id = r1_model_rate_law1.gen_id()
        r1_rate_law_expression2, error = wc_lang.RateLawExpression.deserialize('k_cat_r1_backward_enzyme3 * enzyme3[c]', {
            wc_lang.Parameter: {'k_cat_r1_backward_enzyme3': model.parameters.create(id='k_cat_r1_backward_enzyme3')},
            wc_lang.Species: {'enzyme3[c]': model.species.get_one(id='enzyme3[c]')},
            })
        assert error is None, str(error)
        r1_model_rate_law2 = model.rate_laws.create(
            expression=r1_rate_law_expression2,
            reaction=r1,
            direction=wc_lang.RateLawDirection['backward'])
        r1_model_rate_law2.id = r1_model_rate_law2.gen_id()

        r2 = submodel.reactions.create(id='r2', reversible=True, model=model)           
        r2.participants.append(model.species.get_one(id='m1[c]').species_coefficients.get_or_create(coefficient=-1))
        r2.participants.append(model.species.get_one(id='m2[c]').species_coefficients.get_or_create(coefficient=1))
        r2_rate_law_expression1, error = wc_lang.RateLawExpression.deserialize('k_cat_r2_forward_enzyme1 * enzyme1[c] + k_cat_r2_forward_enzyme2 * enzyme2[c]', {
            wc_lang.Parameter: {'k_cat_r2_forward_enzyme1': model.parameters.create(id='k_cat_r2_forward_enzyme1', value=100.),
                                'k_cat_r2_forward_enzyme2': model.parameters.create(id='k_cat_r2_forward_enzyme2')},
            wc_lang.Species: {'enzyme1[c]': model.species.get_one(id='enzyme1[c]'),
                              'enzyme2[c]': model.species.get_one(id='enzyme2[c]')},
            })
        assert error is None, str(error)
        r2_model_rate_law1 = model.rate_laws.create(
            expression=r2_rate_law_expression1,
            reaction=r2,
            direction=wc_lang.RateLawDirection['forward'])
        r2_model_rate_law1.id = r2_model_rate_law1.gen_id()
        r2_rate_law_expression2, error = wc_lang.RateLawExpression.deserialize('k_cat_r2_backward_enzyme1 * enzyme1[c] + k_cat_r2_backward_enzyme2 * enzyme2[c]', {
            wc_lang.Parameter: {'k_cat_r2_backward_enzyme1': model.parameters.create(id='k_cat_r2_backward_enzyme1'),
                                'k_cat_r2_backward_enzyme2': model.parameters.create(id='k_cat_r2_backward_enzyme2')},
            wc_lang.Species: {'enzyme1[c]': model.species.get_one(id='enzyme1[c]'),
                              'enzyme2[c]': model.species.get_one(id='enzyme2[c]')},
            })
        assert error is None, str(error)
        r2_model_rate_law2 = model.rate_laws.create(
            expression=r2_rate_law_expression2,
            reaction=r2,
            direction=wc_lang.RateLawDirection['backward'])
        r2_model_rate_law2.id = r2_model_rate_law2.gen_id()

        r3 = submodel.reactions.create(id='r3', reversible=False, model=model)           
        r3.participants.append(model.species.get_one(id='m1[c]').species_coefficients.get_or_create(coefficient=-2))
        r3.participants.append(model.species.get_one(id='m3[c]').species_coefficients.get_or_create(coefficient=1))
        r3_rate_law_expression, error = wc_lang.RateLawExpression.deserialize('k_cat_r3_forward_enzyme2 * enzyme2[c]', {
            wc_lang.Parameter: {'k_cat_r3_forward_enzyme2': model.parameters.create(id='k_cat_r3_forward_enzyme2', value=200.)},
            wc_lang.Species: {'enzyme2[c]': model.species.get_one(id='enzyme2[c]')},
            })
        assert error is None, str(error)
        r3_model_rate_law = model.rate_laws.create(
            expression=r3_rate_law_expression,
            reaction=r3,
            direction=wc_lang.RateLawDirection['forward'])
        r3_model_rate_law.id = r3_model_rate_law.gen_id()

        r4 = submodel.reactions.create(id='r4', reversible=False, model=model)           
        r4.participants.append(model.species.get_one(id='m2[c]').species_coefficients.get_or_create(coefficient=-2))
        r4.participants.append(model.species.get_one(id='m3[c]').species_coefficients.get_or_create(coefficient=1))

        biomass_rxn = submodel.reactions.create(id='biomass_reaction', reversible=False, model=model)           
        biomass_rxn.participants.append(model.species.get_one(id='m3[c]').species_coefficients.get_or_create(coefficient=-1))

    def test_gen_reactions(self):

        # Create KB content
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()

        mito = cell.compartments.create(id='m')
        cytoplasm = cell.compartments.create(id='c')

        trans1 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans1')
        trans1_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=trans1, 
            value='10000.0', value_type=wc_ontology['WC:float'])
        trans1_species = wc_kb.core.Species(species_type=trans1, compartment=cytoplasm)
        prot1 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot1', transcript=trans1)
        prot1_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot1, 
            value='20000.0', value_type=wc_ontology['WC:float'])

        trans2 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans2')
        trans2_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=trans2, 
            value='10000.0', value_type=wc_ontology['WC:float'])
        trans2_species = wc_kb.core.Species(species_type=trans2, compartment=cytoplasm)
        prot2 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot2', transcript=trans2)
        prot2_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot2, 
            value='20000.0', value_type=wc_ontology['WC:float'])

        trans3 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans3')
        trans3_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=trans3, 
            value='10000.0', value_type=wc_ontology['WC:float'])
        trans3_species = wc_kb.core.Species(species_type=trans3, compartment=mito)
        prot3 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot3', transcript=trans3)
        prot3_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot3, 
            value='20000.0', value_type=wc_ontology['WC:float'])

        trans4 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans4')
        trans4_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=trans4, 
            value='10000.0', value_type=wc_ontology['WC:float'])
        trans4_species = wc_kb.core.Species(species_type=trans4, compartment=mito)
       
        # Create initial model content
        model = wc_lang.Model()
       
        model.parameters.create(id='mean_doubling_time', value=20000)

        for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType):
            model_species_type = model.species_types.create(id=i.id, type=wc_ontology['WC:RNA'])
            model_compartment = model.compartments.get_or_create(id=i.species[0].compartment.id)
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=model_species, mean=10.)
            conc_model.id = conc_model.gen_id()
            if i.protein:
                model_species_type = model.species_types.create(id=i.protein.id, type=wc_ontology['WC:protein'])                
                model_compartment = model.compartments.get_or_create(id=i.species[0].compartment.id)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, mean=2.)
                conc_model.id = conc_model.gen_id() 
                if i.protein.id == 'prot1':
                    model_compartment = model.compartments.get_or_create(id='n')
                    model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                    model_species.id = model_species.gen_id()
                    conc_model = model.distribution_init_concentrations.create(species=model_species, mean=2.)
                    conc_model.id = conc_model.gen_id()
		            
        metabolic_participants = ['atp', 'ctp', 'gtp', 'utp', 'ppi', 'amp', 'cmp', 'rec', 'pool',
            'gmp', 'ump', 'h2o', 'h', 'adp', 'pi', 'gdp', 'ala_L', 'met_L', 'g6p', 'chsterol']
        for i in metabolic_participants:
            model_species_type = model.species_types.create(id=i, type=wc_ontology['WC:metabolite'])
            for j in ['n', 'm', 'c', 'l']:
                model_compartment = model.compartments.get_or_create(id=j)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
        conc_model = model.distribution_init_concentrations.create(species=model.species.get_one(id='pool[c]'), mean=25.)
        conc_model.id = conc_model.gen_id()        
	            
        others = ['polr2', 'ribosome', 'polr_bound_non_specific_species', 
            'polr_binding_site_species', 'polr_bound_species', 'polr_non_specific_binding_site_species', 
            'ribo_binding_site_species', 'ribo_bound_species']
        for i in others:
            model_species_type = model.species_types.create(id=i, type=wc_ontology['WC:pseudo_species'])
            for j in ['n', 'm', 'c']:
                model_compartment = model.compartments.get_or_create(id=j)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()                        	

        # Create transcription submodel
        transcription_submodel = model.submodels.create(id='transcription')        
        for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType):
            transcription_compartment = 'n' if i.species[0].compartment.id=='c' else 'm' 
            translation_compartment = 'c' if i.species[0].compartment.id=='c' else 'm'
            # Initiation
            init_reaction = model.reactions.create(id='transcription_initiation_' + i.id, submodel=transcription_submodel)
            init_reaction.participants.append(model.species.get_one(
                id='polr_bound_non_specific_species[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(model.species.get_one(
                id='polr_binding_site_species[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(model.species.get_one(
                id='polr_bound_species[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=1))
            init_reaction.participants.append(model.species.get_one(
                id='polr_non_specific_binding_site_species[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=1))
            if i.id == 'trans1':
                init_reaction.participants.append(model.species.get_one(
                    id='atp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                    coefficient=-2))
                init_reaction.participants.append(model.species.get_one(
                    id='h2o[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                    coefficient=-2))
                init_reaction.participants.append(model.species.get_one(
                    id='adp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                    coefficient=2))
                init_reaction.participants.append(model.species.get_one(
                    id='pi[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                    coefficient=2))
                init_reaction.participants.append(model.species.get_one(
                    id='h[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                    coefficient=2))

            # Elongation
            reaction = model.reactions.get_or_create(id='transcription_elongation_' + i.id, submodel=transcription_submodel)
            reaction.participants.append(model.species.get_one(
                id='polr_bound_species[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=-1))
            reaction.participants.append(model.species.get_one(
                id='atp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=-2))
            reaction.participants.append(model.species.get_one(
                id='ctp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=-2))
            reaction.participants.append(model.species.get_one(
                id='gtp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=-2))
            reaction.participants.append(model.species.get_one(
                id='utp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=-2))
            reaction.participants.append(model.species.get_one(
                id='h2o[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=-3))
            reaction.participants.append(model.species.get_one(
                id='{}[{}]'.format(i.id, translation_compartment)).species_coefficients.get_or_create(
                coefficient=1))
            reaction.participants.append(model.species.get_one(
                id='ppi[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=7))
            reaction.participants.append(model.species.get_one(
                id='amp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=2))
            reaction.participants.append(model.species.get_one(
                id='cmp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=2))
            reaction.participants.append(model.species.get_one(
                id='gmp[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=2))
            reaction.participants.append(model.species.get_one(
                id='ump[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=2))
            reaction.participants.append(model.species.get_one(
                id='h[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=3))
            reaction.participants.append(model.species.get_one(
                id='polr2[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=1))
            reaction.participants.append(model.species.get_one(
                id='polr_binding_site_species[{}]'.format(transcription_compartment)).species_coefficients.get_or_create(
                coefficient=1))
            if i.protein:
                reaction.participants.append(model.species.get_one(
                    id='ribo_binding_site_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(
                    coefficient=1))         

        # Create RNA degradation submodel
        rna_deg_submodel = model.submodels.create(id='rna_degradation')
        for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType):
            reaction = model.reactions.get_or_create(id='degradation_' + i.id, submodel=rna_deg_submodel)
            reaction.participants.append(model.species.get_one(
                id='{}[{}]'.format(i.id, i.species[0].compartment.id)).species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.append(model.species.get_one(
                id='h2o[{}]'.format(i.species[0].compartment.id)).species_coefficients.get_or_create(coefficient=-3))
            if i.protein:            
                reaction.participants.append(model.species.get_one(
                    id='ribo_binding_site_species[{}]'.format(i.species[0].compartment.id)).species_coefficients.get_or_create(
                    coefficient=-1))            
            reaction.participants.append(model.species.get_one(
                id='amp[{}]'.format(i.species[0].compartment.id)).species_coefficients.get_or_create(coefficient=1))
            reaction.participants.append(model.species.get_one(
                id='cmp[{}]'.format(i.species[0].compartment.id)).species_coefficients.get_or_create(coefficient=1))
            reaction.participants.append(model.species.get_one(
                id='gmp[{}]'.format(i.species[0].compartment.id)).species_coefficients.get_or_create(coefficient=1))
            reaction.participants.append(model.species.get_one(
                id='ump[{}]'.format(i.species[0].compartment.id)).species_coefficients.get_or_create(coefficient=1))
            reaction.participants.append(model.species.get_one(
                id='h[{}]'.format(i.species[0].compartment.id)).species_coefficients.get_or_create(coefficient=3))

        # Create translation and translocation submodel
        translation_submodel = model.submodels.create(id='translation_translocation')
        for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType):
            translation_compartment = 'c' if i.species[0].compartment.id=='c' else 'm'
            if i.protein:
                # Initiation
                init_reaction = model.reactions.create(id='translation_initiation_' + i.id, submodel=translation_submodel)
                init_reaction.participants.append(model.species.get_one(
                    id='ribosome[{}]'.format(translation_compartment)).species_coefficients.get_or_create(
                    coefficient=-1))
                init_reaction.participants.append(model.species.get_one(
                    id='ribo_binding_site_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(
                    coefficient=-1))
                init_reaction.participants.append(model.species.get_one(
                    id='met_L[{}]'.format(translation_compartment)).species_coefficients.get_or_create(
                    coefficient=-1))
                init_reaction.participants.append(model.species.get_one(
                    id='h2o[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-5))
                init_reaction.participants.append(model.species.get_one(
                    id='atp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-2))
                init_reaction.participants.append(model.species.get_one(
                    id='gtp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-2))
                init_reaction.participants.append(model.species.get_one(
                    id='ribo_bound_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
                init_reaction.participants.append(model.species.get_one(
                    id='h[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=5))
                init_reaction.participants.append(model.species.get_one(
                    id='amp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
                init_reaction.participants.append(model.species.get_one(
                    id='adp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
                init_reaction.participants.append(model.species.get_one(
                    id='gdp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=2))
                init_reaction.participants.append(model.species.get_one(
                    id='pi[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=5))

                # Elongation
                el_reaction = model.reactions.get_or_create(id='translation_elongation_' + i.id, submodel=translation_submodel)
                el_reaction.participants.append(model.species.get_one(
                    id='ribo_bound_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-1))
                el_reaction.participants.append(model.species.get_one(
                    id='gtp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-3))
                el_reaction.participants.append(model.species.get_one(
                    id='atp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-2))
                el_reaction.participants.append(model.species.get_one(
                    id='h2o[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-5))
                el_reaction.participants.append(model.species.get_one(
                    id='ala_L[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-3))            
                el_reaction.participants.append(model.species.get_one(
                    id='ribosome[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
                el_reaction.participants.append(model.species.get_one(
                    id='ribo_binding_site_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
                el_reaction.participants.append(model.species.get_one(
                    id='{}[{}]'.format(i.protein.id, translation_compartment)).species_coefficients.get_or_create(coefficient=1))
                el_reaction.participants.append(model.species.get_one(
                    id='amp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=2))
                el_reaction.participants.append(model.species.get_one(
                    id='gdp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=3))
                el_reaction.participants.append(model.species.get_one(
                    id='pi[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=7))
                el_reaction.participants.append(model.species.get_one(
                    id='h[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=3))
        
                if i.protein.id == 'prot1':
                    # Translocation
                    trans_reaction = model.reactions.get_or_create(id='translocation_prot1_c_to_n', submodel=translation_submodel)
                    trans_reaction.participants.append(model.species.get_one(id='prot1[c]').species_coefficients.get_or_create(
                        coefficient=-1))
                    trans_reaction.participants.append(model.species.get_one(id='gtp[n]').species_coefficients.get_or_create(coefficient=-1))
                    trans_reaction.participants.append(model.species.get_one(id='h2o[n]').species_coefficients.get_or_create(coefficient=-1))
                    trans_reaction.participants.append(model.species.get_one(id='prot1[n]').species_coefficients.get_or_create(coefficient=1))
                    trans_reaction.participants.append(model.species.get_one(id='gdp[n]').species_coefficients.get_or_create(coefficient=1))
                    trans_reaction.participants.append(model.species.get_one(id='pi[n]').species_coefficients.get_or_create(coefficient=1))
                    trans_reaction.participants.append(model.species.get_one(id='h[n]').species_coefficients.get_or_create(coefficient=1))

        # Create protein degradation submodel
        prot_deg_submodel = model.submodels.create(id='protein_degradation')
        degradation_comp = model.compartments.get_one(id='l')
        for protein_model in model.species_types.get(type=wc_ontology['WC:protein']):
            for protein_sp in protein_model.species:
                model_rxn = model.reactions.create(id='{}_{}_degradation'.format(protein_model.id, protein_sp.compartment.id), submodel=prot_deg_submodel)
                model_rxn.participants.add(protein_sp.species_coefficients.get_or_create(coefficient=-1))
                model_rxn.participants.add(model.species.get_one(id='ala_L[l]').species_coefficients.get_or_create(coefficient=3))        
                model_rxn.participants.add(model.species.get_one(id='met_L[l]').species_coefficients.get_or_create(coefficient=1))
                model_rxn.participants.add(model.species.get_one(id='h2o[l]').species_coefficients.get_or_create(coefficient=-3))
                
        gen = metabolism.MetabolismSubmodelGenerator(kb, model, options={
            'recycled_metabolites': {'rec[m]': 100},
            'carbohydrate_components': {'g6p[c]': 1000},
            'lipid_components': {'chsterol[c]': 500},
            'amino_acid_ids': ['ala_L', 'met_L'],
            })
        gen.clean_and_validate_options()
        gen.gen_reactions()
        
        self.assertEqual(model.reactions.get_one(id='biomass_reaction').reversible, False)
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='biomass_reaction').participants},
            {'pool[c]': 25, 'rec[m]': -100, 'g6p[c]': 1000, 'chsterol[c]': 500, 
            'atp[n]': 180, 'ctp[n]': 120, 'gtp[n]': 124, 'utp[n]': 120, 'ppi[n]': -420, 'amp[n]': -120, 'cmp[n]': -120, 'gmp[n]': -120, 
            'ump[n]': -120, 'h2o[n]': 244, 'h[n]': -244, 'adp[n]': -60, 'pi[n]': -64, 'gdp[n]': -4, 
            'atp[m]': 136, 'ctp[m]': 120, 'gtp[m]': 140, 'utp[m]': 120, 'ppi[m]': -420, 'amp[m]': -172, 'cmp[m]': -160, 'gmp[m]': -160, 
            'ump[m]': -160, 'h2o[m]': 340, 'h[m]': -332, 'adp[m]': -4, 'pi[m]': -48, 'gdp[m]': -20, 'ala_L[m]': 12, 'met_L[m]': 4,
            'atp[c]': 48, 'gtp[c]': 60, 'amp[c]': -76, 'cmp[c]': -40, 'gmp[c]': -40, 'ump[c]': -40, 'h2o[c]': 240, 'h[c]': -216, 
            'adp[c]': -12, 'pi[c]': -144, 'gdp[c]': -60, 'ala_L[c]': 36, 'met_L[c]': 12, 'h2o[l]': 24, 'ala_L[l]': -24, 'met_L[l]': -8,
            })

    def test_input_atp_production(self):
        # Create KB content
        kb = wc_kb.KnowledgeBase()
        cell = kb.cell = wc_kb.Cell()

        cytoplasm = cell.compartments.create(id='c')

        trans1 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans1')
        trans1_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=trans1, 
            value='10000.0', value_type=wc_ontology['WC:float'])
        trans1_species = wc_kb.core.Species(species_type=trans1, compartment=cytoplasm)
        prot1 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot1', transcript=trans1)
        prot1_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot1, 
            value='20000.0', value_type=wc_ontology['WC:float'])

        # Create initial model content
        model = wc_lang.Model()
       
        model.parameters.create(id='mean_doubling_time', value=20000)
        
        model_species_type = model.species_types.create(id=trans1.id, type=wc_ontology['WC:RNA'])
        model_compartment = model.compartments.get_or_create(id=trans1.species[0].compartment.id)
        model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
        model_species.id = model_species.gen_id()
        conc_model = model.distribution_init_concentrations.create(species=model_species, mean=10.)
        conc_model.id = conc_model.gen_id()
        
        model_species_type = model.species_types.create(id=trans1.protein.id, type=wc_ontology['WC:protein'])                
        model_compartment = model.compartments.get_or_create(id=trans1.species[0].compartment.id)
        model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
        model_species.id = model_species.gen_id()
        conc_model = model.distribution_init_concentrations.create(species=model_species, mean=2.)
        conc_model.id = conc_model.gen_id() 
                            
        metabolic_participants = ['atp', 'ctp', 'gtp', 'utp', 'ppi', 'amp', 'cmp', 'rec', 'pool',
            'gmp', 'ump', 'h2o', 'h', 'adp', 'pi', 'gdp', 'ala_L', 'met_L', 'g6p', 'chsterol']
        for i in metabolic_participants:
            model_species_type = model.species_types.create(id=i, type=wc_ontology['WC:metabolite'])            
            model_compartment = model.compartments.get_or_create(id='c')
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
        conc_model = model.distribution_init_concentrations.create(species=model.species.get_one(id='pool[c]'), mean=25.)
        conc_model.id = conc_model.gen_id()        
                
        others = ['ribosome', 'ribo_binding_site_species', 'ribo_bound_species']
        for i in others:
            model_species_type = model.species_types.create(id=i, type=wc_ontology['WC:pseudo_species'])
            model_compartment = model.compartments.get_or_create(id='c')
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()                           

        # Create translation and translocation submodel
        translation_submodel = model.submodels.create(id='translation_translocation')        
        translation_compartment = 'c'        
        # Initiation
        init_reaction = model.reactions.create(id='translation_initiation_' + trans1.id, submodel=translation_submodel)
        init_reaction.participants.append(model.species.get_one(
            id='ribosome[{}]'.format(translation_compartment)).species_coefficients.get_or_create(
            coefficient=-1))
        init_reaction.participants.append(model.species.get_one(
            id='ribo_binding_site_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(
            coefficient=-1))
        init_reaction.participants.append(model.species.get_one(
            id='met_L[{}]'.format(translation_compartment)).species_coefficients.get_or_create(
            coefficient=-1))
        init_reaction.participants.append(model.species.get_one(
            id='h2o[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-5))
        init_reaction.participants.append(model.species.get_one(
            id='atp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-2))
        init_reaction.participants.append(model.species.get_one(
            id='gtp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-2))
        init_reaction.participants.append(model.species.get_one(
            id='ribo_bound_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
        init_reaction.participants.append(model.species.get_one(
            id='h[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=5))
        init_reaction.participants.append(model.species.get_one(
            id='amp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
        init_reaction.participants.append(model.species.get_one(
            id='adp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
        init_reaction.participants.append(model.species.get_one(
            id='gdp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=2))
        init_reaction.participants.append(model.species.get_one(
            id='pi[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=5))

        # Elongation
        el_reaction = model.reactions.get_or_create(id='translation_elongation_' + trans1.id, submodel=translation_submodel)
        el_reaction.participants.append(model.species.get_one(
            id='ribo_bound_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-1))
        el_reaction.participants.append(model.species.get_one(
            id='gtp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-3))
        el_reaction.participants.append(model.species.get_one(
            id='atp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-2))
        el_reaction.participants.append(model.species.get_one(
            id='h2o[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-5))
        el_reaction.participants.append(model.species.get_one(
            id='ala_L[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=-3))            
        el_reaction.participants.append(model.species.get_one(
            id='ribosome[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
        el_reaction.participants.append(model.species.get_one(
            id='ribo_binding_site_species[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=1))
        el_reaction.participants.append(model.species.get_one(
            id='{}[{}]'.format(prot1.id, translation_compartment)).species_coefficients.get_or_create(coefficient=1))
        el_reaction.participants.append(model.species.get_one(
            id='amp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=2))
        el_reaction.participants.append(model.species.get_one(
            id='gdp[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=3))
        el_reaction.participants.append(model.species.get_one(
            id='pi[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=7))
        el_reaction.participants.append(model.species.get_one(
            id='h[{}]'.format(translation_compartment)).species_coefficients.get_or_create(coefficient=3))

        gen = metabolism.MetabolismSubmodelGenerator(kb, model, options={
            'amino_acid_ids': ['ala_L', 'met_L'],
            'atp_production': 2000,
            })
        gen.clean_and_validate_options()
        gen.gen_reactions()
        
        self.assertEqual([i.coefficient for i in model.reactions.get_one(id='biomass_reaction').participants if i.species.id=='atp[c]'], [2000])

    def test_calibrate_submodel(self):
        
        model = self.model
        gen = self.gen
        gen.clean_and_validate_options()

        model.species.get_one(id='enzyme3[c]').distribution_init_concentration.mean = 0.01
        model.parameters.get_one(id='k_cat_r1_forward_enzyme3').value = 100.
        model.parameters.get_one(id='k_cat_r1_backward_enzyme3').value = 100.
        model.parameters.get_one(id='k_cat_r2_forward_enzyme2').value = 200.
        model.parameters.get_one(id='k_cat_r2_backward_enzyme1').value = 100.
        model.parameters.get_one(id='k_cat_r2_backward_enzyme2').value = 100.
        model.parameters.get_one(id='k_cat_r3_forward_enzyme2').value = 300.
        
        r4 = model.reactions.get_one(id='r4')
        r4_rate_law_expression, error = wc_lang.RateLawExpression.deserialize('k_cat_r4_forward_enzyme1 * enzyme1[c]', {
            wc_lang.Parameter: {'k_cat_r4_forward_enzyme1': model.parameters.create(id='k_cat_r4_forward_enzyme1', value=600.)},
            wc_lang.Species: {'enzyme1[c]': model.species.get_one(id='enzyme1[c]')},
            })
        assert error is None, str(error)
        r4_model_rate_law = model.rate_laws.create(
            expression=r4_rate_law_expression,
            reaction=r4,
            direction=wc_lang.RateLawDirection['forward'])
        r4_model_rate_law.id = r4_model_rate_law.gen_id()

        gen.options['scale_factor'] = 1e2
        gen.options['coef_scale_factor'] = 10
        gen.options['media_fluxes'] = {'ex_m1': (10., 12.), 'ex_m2': (10., 12.)}
        gen.options['exchange_reactions'] = ['ex_m1', 'ex_m2', 'ex_m3']
        gen.calibrate_submodel()

        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward_enzyme3').value, 300.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward_enzyme3').comments, 'Measured value adjusted to relax bound')
        self.assertEqual(model.parameters.get_one(id='k_cat_r4_forward_enzyme1').value, 600.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r4_forward_enzyme1').comments, '')

    def test_relax_bounds(self):
        
        gen = self.gen
        gen.clean_and_validate_options()
        
        gen._reaction_bounds = {'ex_m1': (10, 10), 'ex_m2': (10, 10), 'ex_m3': (0, 0), 'r1': (0, 1), 'r2': (-2, 3), 'r3': (0, 2), 'r4': (0, None), 'biomass_reaction': (0, None)}        
        lower_bound_adjustable = ['r2']
        upper_bound_adjustable = ['r1', 'r2', 'r3']
        target = {'biomass_reaction': 10}

        alpha_lower, alpha_upper = gen.relax_bounds(target, lower_bound_adjustable, upper_bound_adjustable)

        self.assertEqual(alpha_lower, {})
        self.assertEqual(alpha_upper, {'r3': 1})
    
    def test_flux_variability_analysis(self):
        
        gen = self.gen
        gen.clean_and_validate_options()

        bounds = {'ex_m1': (10, 10), 'ex_m2': (10, 10), 'ex_m3': (0, 0), 'r1': (0, 1), 'r2': (-2, 3), 'r3': (0, 3), 'r4': (0, None), 'biomass_reaction': (0, None)}
        submodel = gen.submodel
        conv_model = conv_opt.Model(name='test_model')
        conv_variables = {}
        conv_metabolite_matrices = collections.defaultdict(list)
        for reaction in submodel.reactions:
            conv_variables[reaction.id] = conv_opt.Variable(
                name=reaction.id, type=conv_opt.VariableType.continuous,
                lower_bound=bounds[reaction.id][0], 
                upper_bound=bounds[reaction.id][1])
            conv_model.variables.append(conv_variables[reaction.id])
            for part in reaction.participants:
                if reaction.id == 'biomass_reaction':
                    conv_metabolite_matrices[part.species.id].append(
                        conv_opt.LinearTerm(conv_variables[reaction.id], 
                            part.coefficient))
                else:
                    conv_metabolite_matrices[part.species.id].append(
                        conv_opt.LinearTerm(conv_variables[reaction.id], 
                            part.coefficient))  

        for met_id, expression in conv_metabolite_matrices.items():
            conv_model.constraints.append(conv_opt.Constraint(expression, name=met_id, 
                upper_bound=0.0, lower_bound=0.0))                      

        conv_model.objective_terms = [conv_opt.LinearTerm(conv_variables['biomass_reaction'], 1.),]
        
        flux_range = gen.flux_variability_analysis(conv_model)
        self.assertEqual(flux_range, {'ex_m1': (10., 10.), 'ex_m2': (10., 10.), 'ex_m3': (0., 0.), 'r1': (1., 1.), 'r2': (3., 3.), 'r3': (3., 3.), 'r4': (6., 6.), 'biomass_reaction': (10., 10.)})

        flux_range = gen.flux_variability_analysis(conv_model, fraction_of_objective=0.8, target_reactions=['r1', 'r2', 'biomass_reaction'])
        self.assertEqual(flux_range, {'r1': (0., 1.), 'r2': (3., 3.), 'biomass_reaction': (8., 8.)})

    def test_impute_kinetic_constant(self): 
        
        model = self.model
        gen = self.gen
        gen.clean_and_validate_options()

        bound_values = {'ex_m1': (10, 10), 'ex_m2': (10, 10), 'ex_m3': (0, 0), 'r1': (0, 1), 'r2': (-2, 3), 'r3': (0, 3), 'r4': (0, None), 'biomass_reaction': (0, None)}
        gen.impute_kinetic_constant(bound_values)

        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward_enzyme3').value, 150.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward_enzyme3').comments, 'Value imputed as the median of measured k_cat values')
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_backward_enzyme3').value, 150.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_backward_enzyme3').comments, 'Value imputed as the median of measured k_cat values')
        self.assertEqual(model.parameters.get_one(id='k_cat_r2_forward_enzyme1').value, 100.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r2_forward_enzyme1').comments, '') 
        self.assertEqual(model.parameters.get_one(id='k_cat_r2_forward_enzyme2').value, 200.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r2_forward_enzyme2').comments, 'Value imputed based on FVA bound value')
        self.assertEqual(model.parameters.get_one(id='k_cat_r2_backward_enzyme1').value, 100.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r2_backward_enzyme1').comments, 'Value imputed based on FVA bound value')
        self.assertEqual(model.parameters.get_one(id='k_cat_r2_backward_enzyme2').value, 100.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r2_backward_enzyme2').comments, 'Value imputed based on FVA bound value') 
        self.assertEqual(model.parameters.get_one(id='k_cat_r3_forward_enzyme2').value, 300.) 
        self.assertEqual(model.parameters.get_one(id='k_cat_r3_forward_enzyme2').comments, 'Measured value adjusted to relax bound')      
