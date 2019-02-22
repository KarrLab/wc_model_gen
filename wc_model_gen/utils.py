""" Utility methods for generating submodels

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-23
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.ontology import wcm_ontology
from wc_utils.util.units import unit_registry
import math
import wc_lang


def calculate_average_synthesis_rate(mean_concentration, half_life, mean_doubling_time):
	""" Calculate the average synthesis rate of a species over a cell cycle

	    Args:
	    	mean_concentration (:obj:`float`): species mean concentration
	    	half_life (:obj:`float`): species half life
	    	mean_doubling_time (:obj:`float`): mean doubling time of cells

	    Returns:
	    	:obj:`float`: the average synthesis rate of the species
	"""
	ave_synthesis_rate = math.log(2) * (1. / mean_doubling_time + 1. / half_life) * mean_concentration

	return ave_synthesis_rate

def calculate_average_degradation_rate(mean_concentration, half_life):
    """ Calculate the average degradation rate of a species over a cell cycle

        Args:
            mean_concentration (:obj:`float`): species mean concentration
            half_life (:obj:`float`): species half life
            
        Returns:
            :obj:`float`: the average degradation rate of the species
    """
    ave_degradation_rate = math.log(2) / half_life * mean_concentration

    return ave_degradation_rate    

def MM_like_rate_law(Avogadro, molecule_units, reaction, modifiers=None, modifier_reactants=None):
    """ Generate a Michaelis-Menten-like rate law. For a multi-substrate reaction,  
        the substrate term is formulated as the multiplication of a Hill equation
        with a coefficient of 1 for each substrate. For multi-steps reaction where
        each step is catalyzed by a different enzyme, the enzyme term is formulated
        as the multiplication of all the enzyme concentrations. 

        Example:

        	Rate = k_cat * [E1] * [E2] * [S1]/(Km_S1 + [S1]) * [S2]/(Km_S2 + [S2])

        	where
        	    k_cat: catalytic constant
            	[En]: concentration of nth enzyme (modifier)
            	[Sn]: concentration of nth substrate
            	Km_Sn: Michaelis-Menten constant for nth substrate   

        Args:
            Avogadro (:obj:`wc_lang.Parameter`): model parameter for Avogadro number
            molecule_units (:obj:`wc_lang.Parameter`): model parameter for molecule units
        	reaction (:obj:`wc_lang.Reaction`): reaction    
        	modifiers (:obj:`list` of :obj:`wc_lang.Observable`): list of observables,
                each of which evaluates to the total concentration of all enzymes that 
                catalyze the same intermediate step in the reaction
            modifier_reactants (:obj:`list` of :obj:`wc_lang.Species`): list of species 
                in modifiers that should be included as reactants in the rate law    
        	
        Returns:
        	:obj:`wc_lang.RateLawExpression`: rate law
        	:obj:`list` of :obj:`wc_lang.Parameter`: list of parameters in the rate law  	
    """
    if modifiers is None:
        modifier_species = []
    else:
        modifier_species = [i for modifier in modifiers for i in modifier.expression.species]

    if modifier_reactants is None:
        additional_reactants = []
    else:
        additional_reactants = modifier_reactants

    parameters = {}
    parameters[Avogadro.id] = Avogadro
    parameters[molecule_units.id] = molecule_units            

    model_k_cat = wc_lang.Parameter(id='k_cat_{}'.format(reaction.id),
                                    type=wcm_ontology['WCM:k_cat'],
                                    units=unit_registry.parse_units('s^-1'))
    parameters[model_k_cat.id] = model_k_cat    

    expression_terms = []
    all_species = {}
    for species in reaction.get_reactants():
                
        if species not in modifier_species or species in additional_reactants:
        
            all_species[species.gen_id()] = species
            
            model_k_m = wc_lang.Parameter(id='K_m_{}_{}'.format(reaction.id, species.species_type.id),
                type=wcm_ontology['WCM:K_m'],                
                units=unit_registry.parse_units('M'))
            parameters[model_k_m.id] = model_k_m
            
            volume = species.compartment.init_density.function_expressions[0].function
            
            expression_terms.append('({} / ({} + {} * {} * {}))'.format(species.gen_id(),
                                                                        species.gen_id(),
                                                                        model_k_m.id, Avogadro.id,
                                                                        volume.id))   

    expression = '{} / {} * {} * {}'.format(
        model_k_cat.id,
        molecule_units.id, 
        ' * '.join([i.id for i in modifiers]), 
        ' * '.join(expression_terms))
    
    rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, {
        wc_lang.Parameter: parameters,
        wc_lang.Species: all_species,
        wc_lang.Observable: {i.id: i for i in modifiers},
        wc_lang.Function: {volume.id: volume},
        })
    assert error is None, str(error)

    return rate_law_expression, list(parameters.values()) 	
