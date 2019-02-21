""" Tests for utility methods

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-02-13
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.ontology import wcm_ontology
from wc_utils.util.units import unit_registry
import wc_model_gen.utils as utils
import math
import scipy.constants
import unittest
import wc_lang


class TestCase(unittest.TestCase):

    def test_calculate_average_synthesis_rate(self):

        test_rate = utils.calculate_average_synthesis_rate(0.5, 300., 36000.)
        self.assertAlmostEqual(test_rate, 0.001164872, places=9)

    def test_MM_like_rate_law(self):

        Avogadro = wc_lang.Parameter(id='Avogadro', value=scipy.constants.Avogadro)
        molecule_units = wc_lang.Parameter(id='molecule_units', value=1.,
            units=unit_registry.parse_units('molecule'))
		
        c = wc_lang.Compartment(id='c', mean_init_volume=0.5)
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

        rate_law, parameters = utils.MM_like_rate_law(
            Avogadro, molecule_units, reaction, modifiers=[modifier1, modifier2], modifier_reactants=[species['s6_c']])

        self.assertEqual(rate_law.expression, 'k_cat_r1 / molecule_units * e1 * e2 * '
            '(s1[c] / (s1[c] + K_m_r1_s1 * Avogadro * volume_c)) * '
            '(s2[c] / (s2[c] + K_m_r1_s2 * Avogadro * volume_c)) * '
            '(s6[c] / (s6[c] + K_m_r1_s6 * Avogadro * volume_c))')
        self.assertEqual(set([i.gen_id() for i in rate_law.species]), set(['s1[c]', 's2[c]', 's6[c]']))
        self.assertEqual(set(rate_law.observables), set([modifier1, modifier2]))
        self.assertEqual(set(rate_law.parameters), set(parameters))        
        self.assertEqual(rate_law.parameters.get_one(id='k_cat_r1').type, wcm_ontology['WCM:k_cat'])
        self.assertEqual(rate_law.parameters.get_one(id='k_cat_r1').units, unit_registry.parse_units('s^-1'))
        self.assertEqual(rate_law.parameters.get_one(id='K_m_r1_s2').type, wcm_ontology['WCM:K_m'])
        self.assertEqual(rate_law.parameters.get_one(id='K_m_r1_s2').value, 0.5/0.5/Avogadro.value)
        self.assertEqual(rate_law.parameters.get_one(id='K_m_r1_s2').units, unit_registry.parse_units('M'))
