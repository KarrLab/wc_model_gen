""" Tests for utility methods

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-02-13
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import wc_model_gen.utils as utils
import math
import scipy.constants
import unittest
import wc_lang


class TestCase(unittest.TestCase):

    def test_calc_avg_syn_rate(self):

        test_rate = utils.calc_avg_syn_rate(0.5, 300., 36000.)
        self.assertAlmostEqual(test_rate, 0.001164872, places=9)

    def test_calc_avg_deg_rate(self):

        test_rate = utils.calc_avg_deg_rate(0.5, 300.)
        self.assertAlmostEqual(test_rate, 0.0011552453009332421, places=16)

    def test_test_metabolite_production(self):

        model = wc_lang.Model()
        submodel = wc_lang.Submodel(model=model, id='metabolism')
        compartment = model.compartments.create(id='c')

        for i in ['o2', 'h2o', 'atp']:
            st = model.species_types.create(id=i)
            species = model.species.create(species_type=st, compartment=compartment)
            species.id = species.gen_id()

        R1 = model.reactions.create(submodel=submodel, id='Ex_o2')
        R1.participants.add(model.species.get_one(
            id='o2[c]').species_coefficients.get_or_create(coefficient=1.0))

        R2 = model.reactions.create(submodel=submodel, id='Ex_h2o')
        R2.participants.add(model.species.get_one(
            id='h2o[c]').species_coefficients.get_or_create(coefficient=1.0))

        R3 = model.reactions.create(submodel=submodel, id='Ex_atp')
        R3.participants.add(model.species.get_one(
            id='atp[c]').species_coefficients.get_or_create(coefficient=-1.0))

        biomass_reaction = model.reactions.create(submodel=submodel, id='biomass_reaction')
        biomass_reaction.participants.add(model.species.get_one(
            id='o2[c]').species_coefficients.get_or_create(coefficient=-1.0))
        biomass_reaction.participants.add(model.species.get_one(
            id='h2o[c]').species_coefficients.get_or_create(coefficient=-1.0))
        biomass_reaction.participants.add(model.species.get_one(
            id='atp[c]').species_coefficients.get_or_create(coefficient=1.0))
        
        submodel.dfba_obj = wc_lang.DfbaObjective(model=model)
        submodel.dfba_obj.id = submodel.dfba_obj.gen_id()
        obj_expression = biomass_reaction.id
        dfba_obj_expression, error = wc_lang.DfbaObjectiveExpression.deserialize(
            obj_expression, {wc_lang.Reaction: {biomass_reaction.id: biomass_reaction}})
        assert error is None, str(error)
        submodel.dfba_obj.expression = dfba_obj_expression

        reaction_bounds = {i.id:(0., 1000.) for i in model.reactions}
        
        unproducibles, unrecyclables = utils.test_metabolite_production(submodel, reaction_bounds, 
            pseudo_reactions=['biomass_reaction'])

        self.assertEqual(unproducibles, [])
        self.assertEqual(unrecyclables, [])

        mock1 = model.species.create(id='mock1')
        mock2 = model.species.create(id='mock2')
        biomass_reaction.participants.add(mock1.species_coefficients.get_or_create(coefficient=1.0))
        biomass_reaction.participants.add(mock2.species_coefficients.get_or_create(coefficient=-1.0))
            
        unproducibles, unrecyclables = utils.test_metabolite_production(submodel, reaction_bounds, 
            pseudo_reactions=['biomass_reaction'])

        self.assertEqual(unproducibles, ['mock1'])
        self.assertEqual(unrecyclables, ['mock2'])

        unproducibles, unrecyclables = utils.test_metabolite_production(submodel, reaction_bounds, 
            test_products=['mock2'], test_reactants=['mock1'])

        self.assertEqual(unproducibles, ['mock2'])
        self.assertEqual(unrecyclables, ['mock1'])

        R4 = model.reactions.create(submodel=submodel, id='Ex_mock1')
        R4.participants.add(mock1.species_coefficients.get_or_create(coefficient=-1.0))

        R5 = model.reactions.create(submodel=submodel, id='Ex_mock2')
        R5.participants.add(mock2.species_coefficients.get_or_create(coefficient=1.0))

        reaction_bounds = {i.id:(0., 1000.) for i in model.reactions}
        
        unproducibles, unrecyclables = utils.test_metabolite_production(submodel, reaction_bounds, 
            test_products=['mock2'], test_reactants=['mock1'])

        self.assertEqual(unproducibles, [])
        self.assertEqual(unrecyclables, [])


    def test_simple_repressor(self):
        model = wc_lang.Model()

        init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], mean=0.5, std=0)
        c = wc_lang.Compartment(id='c', init_volume=init_volume)
        c.init_density = wc_lang.Parameter(id='density_' + c.id, value=1.)

        volume = wc_lang.Function(id='volume_' + c.id)
        volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
        assert error is None, str(error)

        tf_species_type = wc_lang.SpeciesType(id='Repressor')
        tf_species = wc_lang.Species(species_type=tf_species_type, compartment=c)
        tf_species.id = tf_species.gen_id()
        wc_lang.DistributionInitConcentration(species=tf_species, mean=0.5)

        F_rep, species, parameters, functions = utils.simple_repressor(model, 'transcription_rna1', tf_species)

        self.assertEqual(F_rep, '(1 / (1 + Repressor[c] / (Kr_transcription_rna1_Repressor * Avogadro * volume_c)))')
        self.assertEqual(species, {'Repressor[c]': tf_species})
        self.assertEqual(functions, {'volume_c': volume})
        self.assertEqual(set(model.parameters), set(parameters.values()))
        self.assertEqual(sorted(list(parameters.keys())), sorted(['Avogadro', 'Kr_transcription_rna1_Repressor']))
        self.assertEqual(model.parameters.get_one(id='Kr_transcription_rna1_Repressor').type, None)
        self.assertEqual(model.parameters.get_one(id='Kr_transcription_rna1_Repressor').units, unit_registry.parse_units('M'))
     
    def test_simple_activator(self):
        model = wc_lang.Model()

        init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], mean=0.5, std=0)
        c = wc_lang.Compartment(id='c', init_volume=init_volume)
        c.init_density = wc_lang.Parameter(id='density_' + c.id, value=1.)

        volume = wc_lang.Function(id='volume_' + c.id)
        volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
        assert error is None, str(error)

        tf_species_type = wc_lang.SpeciesType(id='Activator')
        tf_species = wc_lang.Species(species_type=tf_species_type, compartment=c)
        tf_species.id = tf_species.gen_id()
        wc_lang.DistributionInitConcentration(species=tf_species, mean=0.5)

        F_act, species, parameters, functions = utils.simple_activator(model, 'transcription_rna1', tf_species)

        self.assertEqual(F_act, '((1 + Activator[c] / (Ka_transcription_rna1_Activator * Avogadro * volume_c) * f_transcription_rna1_Activator) / '
            '(1 + Activator[c] / (Ka_transcription_rna1_Activator * Avogadro * volume_c)))')
        self.assertEqual(species, {'Activator[c]': tf_species})
        self.assertEqual(functions, {'volume_c': volume})
        self.assertEqual(set(model.parameters), set(parameters.values()))
        self.assertEqual(sorted(list(parameters.keys())), sorted(['Avogadro', 'f_transcription_rna1_Activator', 'Ka_transcription_rna1_Activator']))
        self.assertEqual(model.parameters.get_one(id='Ka_transcription_rna1_Activator').type, None)
        self.assertEqual(model.parameters.get_one(id='Ka_transcription_rna1_Activator').units, unit_registry.parse_units('M'))

    def test_gen_michaelis_menten_like_rate_law(self):
        model = wc_lang.Model()

        init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], mean=0.5, std=0)
        c = wc_lang.Compartment(id='c', init_volume=init_volume)
        c.init_density = wc_lang.Parameter(id='density_' + c.id, value=1.)

        volume = wc_lang.Function(id='volume_' + c.id)
        volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
        assert error is None, str(error)

        species_types = {}
        species = {}
        for i in range(1,7):
            Id = 's' + str(i)
            species_types[Id] = wc_lang.SpeciesType(id=Id)
            species[Id + '_c'] = wc_lang.Species(species_type=species_types[Id], compartment=c)
            wc_lang.DistributionInitConcentration(species=species[Id + '_c'], mean=0.5)

        ob_exp1, error = wc_lang.ObservableExpression.deserialize('s4[c] + s5[c]', {
            wc_lang.Species:{species['s4_c'].gen_id(): species['s4_c'],
                            species['s5_c'].gen_id(): species['s5_c']}})
        assert error is None, str(error)
        modifier1 = wc_lang.Observable(id='e1', expression=ob_exp1)

        ob_exp2, error = wc_lang.ObservableExpression.deserialize('2 * s6[c]', {
            wc_lang.Species:{species['s6_c'].gen_id(): species['s6_c']}})
        assert error is None, str(error)
        modifier2 = wc_lang.Observable(id='e2', expression=ob_exp2)

        participant1 = wc_lang.SpeciesCoefficient(species=species['s1_c'], coefficient=-1)
        participant2 = wc_lang.SpeciesCoefficient(species=species['s2_c'], coefficient=-1)
        participant3 = wc_lang.SpeciesCoefficient(species=species['s3_c'], coefficient=1)
        participant4 = wc_lang.SpeciesCoefficient(species=species['s4_c'], coefficient=-1)
        participant5 = wc_lang.SpeciesCoefficient(species=species['s4_c'], coefficient=1)
        participant6 = wc_lang.SpeciesCoefficient(species=species['s6_c'], coefficient=-1)
        participant7 = wc_lang.SpeciesCoefficient(species=species['s6_c'], coefficient=-1)
        participant8 = wc_lang.SpeciesCoefficient(species=species['s6_c'], coefficient=1)
        reaction = wc_lang.Reaction(id='r1', participants=[participant1, participant2, participant3,
            participant4, participant5, participant6, participant7, participant8])

        rate_law, parameters = utils.gen_michaelis_menten_like_rate_law(
            model, reaction, modifiers=[modifier1, modifier2], modifier_reactants=[species['s6_c']])

        self.assertEqual(rate_law.expression, 'k_cat_r1 * e1 * e2 * '
            '(s1[c] / (s1[c] + K_m_r1_s1 * Avogadro * volume_c)) * '
            '(s2[c] / (s2[c] + K_m_r1_s2 * Avogadro * volume_c)) * '
            '(s6[c] / (s6[c] + K_m_r1_s6 * Avogadro * volume_c))')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set(['s1[c]', 's2[c]', 's6[c]']))
        self.assertEqual(set(rate_law.observables), set([modifier1, modifier2]))
        self.assertEqual(set(rate_law.parameters), set(parameters))
        self.assertEqual(rate_law.parameters.get_one(id='k_cat_r1').type, wc_ontology['WC:k_cat'])
        self.assertEqual(rate_law.parameters.get_one(id='k_cat_r1').units, unit_registry.parse_units('s^-1 molecule^-2'))
        self.assertEqual(rate_law.parameters.get_one(id='K_m_r1_s2').type, wc_ontology['WC:K_m'])
        self.assertEqual(rate_law.parameters.get_one(id='K_m_r1_s2').units, unit_registry.parse_units('M'))

        reaction = wc_lang.Reaction(id='r1', participants=[participant1, participant2, participant4, participant8])
        rate_law, parameters = utils.gen_michaelis_menten_like_rate_law(
            model, reaction)
        self.assertEqual(rate_law.expression, 'k_cat_r1 * '
            '(s1[c] / (s1[c] + K_m_r1_s1 * Avogadro * volume_c)) * '
            '(s2[c] / (s2[c] + K_m_r1_s2 * Avogadro * volume_c)) * '
            '(s4[c] / (s4[c] + K_m_r1_s4 * Avogadro * volume_c))')

        reaction = wc_lang.Reaction(id='r1', participants=[participant3])
        rate_law, parameters = utils.gen_michaelis_menten_like_rate_law(
            model, reaction)
        self.assertEqual(rate_law.expression, 'k_cat_r1')

        reaction = wc_lang.Reaction(id='r1', participants=[participant3, participant6])
        rate_law, parameters = utils.gen_michaelis_menten_like_rate_law(
            model, reaction, modifiers=[modifier1, species['s6_c']])
        self.assertEqual(rate_law.expression, 'k_cat_r1 * e1 * s6[c]')

        reaction = wc_lang.Reaction(id='r1', participants=[participant1, participant2, participant4, participant8])
        rate_law, parameters = utils.gen_michaelis_menten_like_rate_law(
            model, reaction, exclude_substrates=[species['s1_c']])
        self.assertEqual(rate_law.expression, 'k_cat_r1 * '
            '(s2[c] / (s2[c] + K_m_r1_s2 * Avogadro * volume_c)) * '
            '(s4[c] / (s4[c] + K_m_r1_s4 * Avogadro * volume_c))')

        with self.assertRaises(TypeError) as ctx:
            rate_law, parameters = utils.gen_michaelis_menten_like_rate_law(
                model, reaction, modifiers=['s6_c'])
        self.assertEqual('The modifiers contain element(s) that is not an observable or a species', str(ctx.exception))  

    def test_gen_michaelis_menten_like_propensity_function(self):
        model = wc_lang.Model()

        init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], mean=0.5, std=0)
        c = wc_lang.Compartment(id='c', init_volume=init_volume)
        c.init_density = wc_lang.Parameter(id='density_' + c.id, value=1.)

        volume = wc_lang.Function(id='volume_' + c.id)
        volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
        assert error is None, str(error)

        species_types = {}
        species = {}
        for i in range(1,7):
            Id = 's' + str(i)
            species_types[Id] = wc_lang.SpeciesType(id=Id)
            model_species = wc_lang.Species(species_type=species_types[Id], compartment=c)
            model_species.id = model_species.gen_id()
            species[Id + '_c'] = model_species 
            wc_lang.DistributionInitConcentration(species=species[Id + '_c'], mean=0.5)

        participant1 = wc_lang.SpeciesCoefficient(species=species['s1_c'], coefficient=-1)
        participant2 = wc_lang.SpeciesCoefficient(species=species['s2_c'], coefficient=-1)
        participant3 = wc_lang.SpeciesCoefficient(species=species['s3_c'], coefficient=-1)
        participant4 = wc_lang.SpeciesCoefficient(species=species['s4_c'], coefficient=-1)
        participant5 = wc_lang.SpeciesCoefficient(species=species['s5_c'], coefficient=1)
        participant6 = wc_lang.SpeciesCoefficient(species=species['s6_c'], coefficient=1)
        
        reaction = wc_lang.Reaction(id='r1', participants=[participant1, participant2, participant3,
            participant4, participant5, participant6])

        with self.assertRaises(ValueError):
            rate_law1, parameters = utils.gen_michaelis_menten_like_propensity_function(
                model, reaction)

        rate_law2, parameters = utils.gen_michaelis_menten_like_propensity_function(
            model, reaction, substrates_as_modifiers=[species['s3_c']])
        self.assertEqual(rate_law2.expression, 'k_cat_r1 * s3[c] * '
            '(s1[c] / (s1[c] + K_m_r1_s1 * Avogadro * volume_c)) * '
            '(s2[c] / (s2[c] + K_m_r1_s2 * Avogadro * volume_c)) * '
            '(s4[c] / (s4[c] + K_m_r1_s4 * Avogadro * volume_c))') 
        self.assertEqual(set([i.gen_id() for i in rate_law2.species]), set(['s1[c]', 's2[c]', 's3[c]', 's4[c]']))
        self.assertEqual(set(rate_law2.parameters), set(parameters))
        self.assertEqual(rate_law2.parameters.get_one(id='k_cat_r1').type, wc_ontology['WC:k_cat'])
        self.assertEqual(rate_law2.parameters.get_one(id='k_cat_r1').units, unit_registry.parse_units('s^-1 molecule^-1'))
        self.assertEqual(rate_law2.parameters.get_one(id='K_m_r1_s2').type, wc_ontology['WC:K_m'])
        self.assertEqual(rate_law2.parameters.get_one(id='K_m_r1_s2').units, unit_registry.parse_units('M'))

        rate_law3, parameters = utils.gen_michaelis_menten_like_propensity_function(
            model, reaction, substrates_as_modifiers=[species['s3_c']], exclude_substrates=[species['s1_c']])
        self.assertEqual(rate_law3.expression, 'k_cat_r1 * s3[c] * '
            '(s2[c] / (s2[c] + K_m_r1_s2 * Avogadro * volume_c)) * '
            '(s4[c] / (s4[c] + K_m_r1_s4 * Avogadro * volume_c))') 
        self.assertEqual(set([i.gen_id() for i in rate_law3.species]), set(['s2[c]', 's3[c]', 's4[c]']))
        self.assertEqual(set(rate_law3.parameters), set(parameters))

    def test_gen_response_functions(self):
        model = wc_lang.Model()
        beta = 2

        init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], mean=0.5, std=0)
        c = wc_lang.Compartment(id='c', name='cytosol', init_volume=init_volume)
        c.init_density = wc_lang.Parameter(id='density_' + c.id, value=1.)

        volume = wc_lang.Function(id='volume_' + c.id)
        volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
        assert error is None, str(error)

        reaction = wc_lang.Reaction()

        species_types = {}
        species = {}
        for i in range(1,5):
            Id = 's' + str(i)
            species_types[Id] = wc_lang.SpeciesType(model=model, id=Id, name='species_type_{}'.format(i))
            model_species = wc_lang.Species(model=model, species_type=species_types[Id], compartment=c)
            model_species.id = model_species.gen_id()
            species[Id + '_c'] = model_species 
            wc_lang.DistributionInitConcentration(species=species[Id + '_c'], mean=0.5)

        factors = [['s1', 'species_type_2'], ['s3'], ['species_type_4']] 
        factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
            model, beta, 'reaction_id', 'reaction_class', c, factors)    

        self.assertEqual(factor_exp, [
            '(reaction_class_factors_c_1 / (reaction_class_factors_c_1 + K_m_reaction_class_reaction_class_factors_c_1 * Avogadro * volume_c))', 
            '(s3[c] / (s3[c] + K_m_reaction_id_s3 * Avogadro * volume_c))', 
            '(s4[c] / (s4[c] + K_m_reaction_id_s4 * Avogadro * volume_c))'])
        self.assertEqual(all_species, {'s1[c]': species['s1_c'], 's2[c]': species['s2_c'], 's3[c]': species['s3_c'], 's4[c]': species['s4_c']})
        self.assertEqual(len(all_parameters), 4)
        self.assertEqual(all_parameters['Avogadro'].value, scipy.constants.Avogadro)
        self.assertEqual(all_parameters['Avogadro'].units, unit_registry.parse_units('molecule mol^-1'))
        self.assertEqual(all_parameters['K_m_reaction_class_reaction_class_factors_c_1'].value, beta * 1. / scipy.constants.Avogadro / c.init_volume.mean)
        self.assertEqual(all_parameters['K_m_reaction_class_reaction_class_factors_c_1'].comments, 'The value was assumed to be 2 times the value of reaction_class_factors_c_1')
        self.assertEqual(all_parameters['K_m_reaction_id_s3'].value, beta * 0.5 / scipy.constants.Avogadro / c.init_volume.mean)
        self.assertEqual(all_parameters['K_m_reaction_id_s4'].type, wc_ontology['WC:K_m'])
        self.assertEqual(all_parameters['K_m_reaction_id_s4'].units, unit_registry.parse_units('M'))
        self.assertEqual(all_parameters['K_m_reaction_id_s4'].comments, 'The value was assumed to be 2 times the concentration of s4 in cytosol')
        self.assertEqual(all_volumes, {'volume_c': volume})
        self.assertEqual(len(all_observables), 1)
        self.assertEqual(len(model.observables), 1)
        self.assertEqual(all_observables['reaction_class_factors_c_1'].name, 'factor for reaction_class in cytosol')
        self.assertEqual(all_observables['reaction_class_factors_c_1'].units, unit_registry.parse_units('molecule'))
        self.assertEqual(all_observables['reaction_class_factors_c_1'].expression.expression, 's1[c] + s2[c]')

        for i in range(5,9):
            Id = 's' + str(i)
            species_types[Id] = wc_lang.SpeciesType(model=model, id=Id, name='species_type_{}'.format(i))
            model_species = wc_lang.Species(model=model, species_type=species_types[Id], compartment=c)
            model_species.id = model_species.gen_id()
            species[Id + '_c'] = model_species 
            wc_lang.DistributionInitConcentration(species=species[Id + '_c'], mean=0.)

        factors = [['s5', 'species_type_6'], ['s7'], ['species_type_8']] 
        factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
            model, beta, 'reaction_id', 'reaction_class', c, factors)

        self.assertEqual(len(model.observables), 2)
        self.assertEqual(all_parameters['K_m_reaction_class_reaction_class_factors_c_2'].value, 1e-05)
        self.assertEqual(all_parameters['K_m_reaction_class_reaction_class_factors_c_2'].comments, 
            'The value was assigned to 1e-05 because the value of reaction_class_factors_c_2 was zero')
        self.assertEqual(all_parameters['K_m_reaction_id_s7'].value, 1e-05)
        self.assertEqual(all_parameters['K_m_reaction_id_s8'].comments, 'The value was assigned to 1e-05 because the concentration of s8 in cytosol was zero')

        factors = [['s5', 'species_type_6']] 
        factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
            model, beta, 'reaction_id2', 'reaction_class2', c, factors)

        self.assertEqual(len(model.observables), 2)
        self.assertEqual(all_parameters['K_m_reaction_class2_reaction_class_factors_c_2'].value, 1e-05)

    def test_gen_mass_action_rate_law(self):
        model = wc_lang.Model()
        c = wc_lang.Compartment(id='c', init_volume=wc_lang.InitVolume(mean=0.5))
        c.init_density = wc_lang.Parameter(id='density_' + c.id, value=1.) 
        kinetic_parameter = wc_lang.Parameter(id='this_parameter', value=1.) 
        volume = wc_lang.Function(id='volume_' + c.id)
        volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
        assert error is None, str(error)
        
        species_types = {}
        species = {}
        for i in range(1,7):
            Id = 's' + str(i)
            species_types[Id] = wc_lang.SpeciesType(id=Id)
            species[Id + '_c'] = wc_lang.Species(species_type=species_types[Id], compartment=c)
            wc_lang.DistributionInitConcentration(species=species[Id + '_c'], mean=0.5)
        Id = 'e'
        species_types[Id] = wc_lang.SpeciesType(id=Id)
        species[Id + '_c'] = wc_lang.Species(species_type=species_types[Id], compartment=c)
        wc_lang.DistributionInitConcentration(species=species[Id + '_c'], mean=0.5)

        # ob_exp1, error = wc_lang.ObservableExpression.deserialize('s4[c] + s5[c]', {
        #     wc_lang.Species:{species['s4_c'].gen_id(): species['s4_c'], 
        #                     species['s5_c'].gen_id(): species['s5_c']}})
        # assert error is None, str(error)
        # modifier1 = wc_lang.Observable(id='e1', expression=ob_exp1)

        # ob_exp2, error = wc_lang.ObservableExpression.deserialize('2 * s6[c]', {
        #     wc_lang.Species:{species['s6_c'].gen_id(): species['s6_c']}})
        # assert error is None, str(error)
        # modifier2 = wc_lang.Observable(id='e2', expression=ob_exp2)
            
        participant1 = wc_lang.SpeciesCoefficient(species=species['s1_c'], coefficient=-1)
        participant2 = wc_lang.SpeciesCoefficient(species=species['s2_c'], coefficient=-1)
        participant3 = wc_lang.SpeciesCoefficient(species=species['s3_c'], coefficient=1)
        participant4 = wc_lang.SpeciesCoefficient(species=species['s4_c'], coefficient=1)
        enzyme_lhs = wc_lang.SpeciesCoefficient(species=species['e_c'], coefficient=-1)
        enzyme_rhs = wc_lang.SpeciesCoefficient(species=species['e_c'], coefficient=1)
        
        reaction = wc_lang.Reaction(id='Assosication', participants=[participant1, participant2, participant3])
        rate_law, parameters = utils.gen_mass_action_rate_law(model, reaction, kinetic_parameter)
        self.assertTrue(rate_law.expression == 'this_parameter * s1[c] * s2[c]' or
            rate_law.expression == 'this_parameter * s2[c] * s1[c]')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set(['s1[c]', 's2[c]']))
        # self.assertEqual(set(rate_law.observables), set([modifier1, modifier2]))
        self.assertEqual(set(rate_law.parameters), set(parameters))        
        # self.assertEqual(rate_law.parameters.get_one(id='k_r1').type, wc_ontology['WC:k_cat'])
        self.assertEqual(rate_law.parameters.get_one(id='this_parameter').units, unit_registry.parse_units('s^-1 * molecule^-1'))
        # self.assertEqual(rate_law.parameters.get_one(id='K_m_r1_s2').type, wc_ontology['WC:K_m'])
        # self.assertEqual(rate_law.parameters.get_one(id='K_m_r1_s2').units, unit_registry.parse_units('M'))

        reaction = wc_lang.Reaction(id='Dissociation', participants=[participant1, participant3, participant4])
        rate_law, parameters = utils.gen_mass_action_rate_law(model, reaction, kinetic_parameter)
        self.assertEqual(rate_law.expression, 'this_parameter * s1[c]')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set(['s1[c]']))     
        self.assertEqual(rate_law.parameters.get_one(id='this_parameter').units, unit_registry.parse_units('s^-1'))

        reaction = wc_lang.Reaction(id='Degradation1', participants=[participant1])
        rate_law, parameters = utils.gen_mass_action_rate_law(model, reaction, kinetic_parameter)
        self.assertEqual(rate_law.expression, 'this_parameter * s1[c]')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set(['s1[c]']))     
        self.assertEqual(rate_law.parameters.get_one(id='this_parameter').units, unit_registry.parse_units('s^-1'))
       
        reaction = wc_lang.Reaction(id='Degradation2', participants=[participant1, enzyme_lhs, enzyme_rhs])
        rate_law, parameters = utils.gen_mass_action_rate_law(model, reaction, kinetic_parameter)
        self.assertTrue(rate_law.expression, 'this_parameter * s1[c] * e[c]' or
            'this_parameter * e[c] * s1[c]')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set(['s1[c]', 'e[c]']))     
        self.assertEqual(rate_law.parameters.get_one(id='this_parameter').units, unit_registry.parse_units('s^-1 * molecule^-1'))

        reaction = wc_lang.Reaction(id='Synthesis1', participants=[participant3])
        rate_law, parameters = utils.gen_mass_action_rate_law(model, reaction, kinetic_parameter)
        self.assertEqual(rate_law.expression, 'this_parameter')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set([]))     
        self.assertEqual(rate_law.parameters.get_one(id='this_parameter').units, unit_registry.parse_units('s^-1 * molecule'))

        reaction = wc_lang.Reaction(id='Synthesis2', participants=[enzyme_lhs, enzyme_rhs, participant3])
        rate_law, parameters = utils.gen_mass_action_rate_law(model, reaction, kinetic_parameter)
        self.assertEqual(rate_law.expression, 'this_parameter * e[c]')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set(['e[c]']))     
        self.assertEqual(rate_law.parameters.get_one(id='this_parameter').units, unit_registry.parse_units('s^-1'))

        reaction = wc_lang.Reaction(id='Conversion', participants=[participant1, enzyme_lhs, enzyme_rhs, participant3]) # Ask Yin Hoon why I can add as many copies of participant2 as I want.
        rate_law, parameters = utils.gen_mass_action_rate_law(model, reaction, kinetic_parameter)
        self.assertTrue(rate_law.expression == 'this_parameter * s1[c] * e[c]' or
            rate_law.expression == 'this_parameter * e[c] * s1[c]')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set(['s1[c]', 'e[c]']))
        self.assertEqual(rate_law.parameters.get_one(id='this_parameter').units, unit_registry.parse_units('s^-1 * molecule^-1'))