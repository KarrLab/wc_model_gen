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

        self.assertEqual(F_rep, '(1 / (1 + Repressor[c] / (Kr_transcription_rna1_Repressor[c] * Avogadro * volume_c))')
        self.assertEqual(species, {'Repressor[c]': tf_species})
        self.assertEqual(functions, {'volume_c': volume})
        self.assertEqual(set(model.parameters), set(parameters.values()))
        self.assertEqual(model.parameters.get_one(id='Kr_transcription_rna1_Repressor[c]').type, None)
        self.assertEqual(model.parameters.get_one(id='Kr_transcription_rna1_Repressor[c]').units, unit_registry.parse_units('M'))
     
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

        self.assertEqual(F_act, '(1 + Activator[c] / (Ka_transcription_rna1_Activator[c] * Avogadro * volume_c) * f_transcription_rna1_Activator[c]) / '
            '(1 + Activator[c] / (Ka_transcription_rna1_Activator[c] * Avogadro * volume_c))')
        self.assertEqual(species, {'Activator[c]': tf_species})
        self.assertEqual(functions, {'volume_c': volume})
        self.assertEqual(set(model.parameters), set(parameters.values()))
        self.assertEqual(model.parameters.get_one(id='Ka_transcription_rna1_Activator[c]').type, None)
        self.assertEqual(model.parameters.get_one(id='Ka_transcription_rna1_Activator[c]').units, unit_registry.parse_units('M'))

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