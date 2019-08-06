""" Utility methods for generating submodels

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-23
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import math
import scipy.constants
import wc_lang


def calc_avg_syn_rate(mean_concentration, half_life, mean_doubling_time):
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


def calc_avg_deg_rate(mean_concentration, half_life):
    """ Calculate the average degradation rate of a species over a cell cycle

        Args:
            mean_concentration (:obj:`float`): species mean concentration
            half_life (:obj:`float`): species half life

        Returns:
            :obj:`float`: the average degradation rate of the species
    """
    ave_degradation_rate = math.log(2) / half_life * mean_concentration

    return ave_degradation_rate


def simple_repressor (model, reaction_id, repressor):
    """ Generate the parameters and string expression of the regulation factor 
        derived in Bintu et al (2005) for the case of a simple repressor

        Args:
            model (:obj:`wc_lang.Model`): model
            reaction_id (:obj:`str`): reaction id   
            repressor (:obj:`wc_lang.Species`): repressor

        Returns:
            :obj:`str`: string expression of the regulation factor
            :obj:`dict` of :obj:`wc_lang.Species`: dict of species in the expression 
                with the species ids as keys and the species objects as values
            :obj:`dict` of :obj:`wc_lang.Parameter`: dict of parameters in the expression 
                with the parameter ids as keys and the parameter objects as values
            :obj:`dict` of :obj:`wc_lang.Function`: dict of functions in the expression 
                with the function ids as keys and the function objects as values                     
    """  
    species = {}
    parameters = {}
    functions = {}

    species[repressor.id] = repressor

    avogadro = model.parameters.get_or_create(
        id='Avogadro',
        type=None,
        value=scipy.constants.Avogadro,
        units=unit_registry.parse_units('molecule mol^-1'))
    parameters[avogadro.id] = avogadro

    Kr = model.parameters.get_or_create(
            id='Kr_{}_{}'.format(reaction_id, repressor.species_type.id),
            type=None,
            units=unit_registry.parse_units('M'))
    parameters[Kr.id] = Kr

    volume = repressor.compartment.init_density.function_expressions[0].function
    functions[volume.id] = volume

    F_rep = '(1 / (1 + {} / ({} * {} * {})))'.format(repressor.id, Kr.id, avogadro.id, volume.id)

    return F_rep, species, parameters, functions


def simple_activator (model, reaction_id, activator):
    """ Generate the parameters and string expression of the regulation factor 
        derived in Bintu et al (2005) for the case of a simple activator

        Args:
            model (:obj:`wc_lang.Model`): model
            reaction_id (:obj:`str`): reaction id
            activator (:obj:`wc_lang.Species`): activator

        Returns:
            :obj:`str`: string expression of the regulation factor
            :obj:`dict` of :obj:`wc_lang.Species`: dict of species in the expression 
                with the species ids as keys and the species objects as values
            :obj:`dict` of :obj:`wc_lang.Parameter`: dict of parameters in the expression 
                with the parameter ids as keys and the parameter objects as values
            :obj:`dict` of :obj:`wc_lang.Function`: dict of functions in the expression 
                with the function ids as keys and the function objects as values            
    """  
    species = {}
    parameters = {}
    functions = {}

    species[activator.id] = activator

    avogadro = model.parameters.get_or_create(
        id='Avogadro',
        type=None,
        value=scipy.constants.Avogadro,
        units=unit_registry.parse_units('molecule mol^-1'))
    parameters[avogadro.id] = avogadro

    Ka = model.parameters.get_or_create(
        id='Ka_{}_{}'.format(reaction_id, activator.species_type.id),
        type=None,
        units=unit_registry.parse_units('M'))
    parameters[Ka.id] = Ka
    
    f = model.parameters.get_or_create(
        id='f_{}_{}'.format(reaction_id, activator.species_type.id),
        type=None)
    parameters[f.id] = f

    volume = activator.compartment.init_density.function_expressions[0].function
    functions[volume.id] = volume
    
    F_act = '((1 + {} / ({} * {} * {}) * {}) / (1 + {} / ({} * {} * {})))'.format(
        activator.id, Ka.id, avogadro.id, volume.id, f.id, activator.id, Ka.id, avogadro.id, volume.id)   

    return F_act, species, parameters, functions


def gen_michaelis_menten_like_rate_law(model, reaction, modifiers=None, modifier_reactants=None):
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
            model (:obj:`wc_lang.Model`): model
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

    avogadro = model.parameters.get_or_create(
        id='Avogadro',
        type=None,
        value=scipy.constants.Avogadro,
        units=unit_registry.parse_units('molecule mol^-1'))
    parameters[avogadro.id] = avogadro

    model_k_cat = model.parameters.get_or_create(id='k_cat_{}'.format(reaction.id),
                                                 type=wc_ontology['WC:k_cat'],
                                                 units=unit_registry.parse_units('s^-1{}'.format(
                                                    (' * molecule^{{-{}}}'.format(len(modifiers))) if modifiers else '')))
    parameters[model_k_cat.id] = model_k_cat

    expression_terms = []
    all_species = {}
    all_volumes = {}
    for species in reaction.get_reactants():

        if species not in modifier_species or species in additional_reactants:

            all_species[species.gen_id()] = species

            model_k_m = model.parameters.create(id='K_m_{}_{}'.format(reaction.id, species.species_type.id),
                                                type=wc_ontology['WC:K_m'],
                                                units=unit_registry.parse_units('M'))
            parameters[model_k_m.id] = model_k_m

            volume = species.compartment.init_density.function_expressions[0].function
            all_volumes[volume.id] = volume

            expression_terms.append('({} / ({} + {} * {} * {}))'.format(species.gen_id(),
                                                                        species.gen_id(),
                                                                        model_k_m.id, avogadro.id,
                                                                        volume.id))

    expression = '{}{}{}'.format(
        model_k_cat.id,
        (' * {}'.format(' * '.join([i.id for i in modifiers]))) if modifiers else '',
        (' * {}'.format(' * '.join(expression_terms))) if expression_terms else '')

    rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, {
        wc_lang.Parameter: parameters,
        wc_lang.Species: all_species,
        wc_lang.Observable: {i.id: i for i in modifiers} if modifiers else {},
        wc_lang.Function: all_volumes,
    })
    assert error is None, str(error)

    return rate_law_expression, list(parameters.values())


def gen_michaelis_menten_like_propensity_function(model, reaction, substrates_as_modifiers=None, exclude_substrates=None):
    """ Generate a Michaelis-Menten-like propensity function. 
        For species that are considered 'substrates', the substrate term is formulated as the 
        multiplication of a Hill equation with a coefficient of 1 for each 'substrate'. 
        For species that are considered 'modifiers', the modifier term is formulated
        as the multiplication of the modifier concentrations. 

        Example:

                Rate = k_cat * [E1] * [E2] * [S1]/(Km_S1 + [S1]) * [S2]/(Km_S2 + [S2])

                where
                    k_cat: catalytic constant
                    [En]: concentration of nth modifier
                    [Sn]: concentration of nth substrate
                    Km_Sn: Michaelis-Menten constant for nth substrate   

        Args:
            model (:obj:`wc_lang.Model`): model
            reaction (:obj:`wc_lang.Reaction`): reaction    
            substrates_as_modifiers (:obj:`list` of :obj:`wc_lang.Species`): list of reactant species 
                that should be considered as modifiers in the rate law
            exclude_substrates (:obj:`list` of :obj:`wc_lang.Species`): list of reactant species 
                that would be excluded from the rate law        

        Returns:
                :obj:`wc_lang.RateLawExpression`: rate law
                :obj:`list` of :obj:`wc_lang.Parameter`: list of parameters in the rate law     
    """
    if substrates_as_modifiers is None:
        raise ValueError('No list has been provided for the input argument substrates_as_modifiers')

    if exclude_substrates:
        excluded_reactants = exclude_substrates
    else:
        excluded_reactants = []        
    
    parameters = {}

    avogadro = model.parameters.get_or_create(
        id='Avogadro',
        type=None,
        value=scipy.constants.Avogadro,
        units=unit_registry.parse_units('molecule mol^-1'))
    parameters[avogadro.id] = avogadro

    model_k_cat = model.parameters.get_or_create(id='k_cat_{}'.format(reaction.id),
                                                 type=wc_ontology['WC:k_cat'],
                                                 units=unit_registry.parse_units('s^-1{}'.format(
                                                    (' * molecule^{{-{}}}'.format(
                                                    len(substrates_as_modifiers))))))
    parameters[model_k_cat.id] = model_k_cat

    expression_terms = []
    all_species = {}
    all_volumes = {}
    for species in reaction.get_reactants():

        all_species[species.gen_id()] = species

        if species not in substrates_as_modifiers and species not in excluded_reactants:            

            model_k_m = model.parameters.create(id='K_m_{}_{}'.format(reaction.id, species.species_type.id),
                                                type=wc_ontology['WC:K_m'],
                                                units=unit_registry.parse_units('M'))
            parameters[model_k_m.id] = model_k_m

            volume = species.compartment.init_density.function_expressions[0].function
            all_volumes[volume.id] = volume

            expression_terms.append('({} / ({} + {} * {} * {}))'.format(species.gen_id(),
                                                                        species.gen_id(),
                                                                        model_k_m.id, avogadro.id,
                                                                        volume.id))

    expression = '{}{}{}'.format(
        model_k_cat.id,
        (' * {}'.format(' * '.join([i.id for i in substrates_as_modifiers]))),
        (' * {}'.format(' * '.join(expression_terms))))

    rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, {
        wc_lang.Parameter: parameters,
        wc_lang.Species: all_species,
        wc_lang.Function: all_volumes,
    })
    assert error is None, str(error)

    return rate_law_expression, list(parameters.values())    


def gen_mass_action_rate_law(model, reaction, model_k, modifiers=None, modifier_reactants=None):
    """ Generate a mass action rate law.

        Example:

            Rate = k * [E1] * [S1]

            where
                k_: rate constant (e.g.: association, dissociation or catalytic constant)
                [En]: concentration of nth enzyme (modifier)
                [Sn]: concentration of nth substrate  

        Args:
            model (:obj:`wc_lang.Model`): model
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
    print('the test reactants are')
    print(reaction.get_reactants())
    print(reaction.get_products())

    if modifiers is None:
        modifier_species = []
        modifier_product = ''
        modifier_ids = None
    else:
        modifier_species = [i for modifier in modifiers for i in modifier.expression.species]
        modifier_product = ' * ' + ' * '.join([i.id for i in modifiers])
        modifier_ids = {i.id: i for i in modifiers}

    if modifier_reactants is None:
        additional_reactants = []
    else:
        additional_reactants = modifier_reactants

    parameters = {}
    
    model_k_unit = len(reaction.get_reactants())-1
    if model_k_unit == 0:
        model_k_unit = 's^-1'
    elif model_k_unit == -1:
        model_k_unit = 's^-1 * molecule'
    else:
        model_k_unit = 's^-1 * molecule^{}'.format(-model_k_unit)
    print('do i have paramers')
    print(len(model.parameters))
    for i in range(len(model.parameters)):
        print(model.parameters[i].id)
    print('i want to get the parameter:')

                                             # type=None,
                                             # units=unit_registry.parse_units(model_k_unit),
                                             # value=model.parameters.get_or_create(id = kinetic_param_name))
    print('i got the parameter:')
    print(model_k.id)
    print(model_k.value)
    # model_k.type = None
    model_k.units = unit_registry.parse_units(model_k_unit)

    parameters[model_k.id] = model_k

    expression_terms = []
    all_species = {}
    all_volumes = {}    
    for species in set(reaction.get_reactants()).union(set(reaction.get_products())):
        if species not in modifier_species or species in additional_reactants:
            all_species[species.gen_id()] = species
            print('added species: {}'.format(species))
            volume = species.compartment.init_density.function_expressions[0].function
            all_volumes[volume.id] = volume
        if species in set(reaction.get_reactants()):
            expression_terms.append(str(species.gen_id()))
    print('the model_k id is:')
    print(model_k.id)
    print('all species is: {}'.format(all_species))
    print('the expression terms are: {}'.format(expression_terms))

    reactant_product = ''
    if expression_terms:
        reactant_product = ' * ' + ' * '.join(expression_terms)
    expression = model_k.id + modifier_product + reactant_product
    print('the expression is')
    print(expression)

    rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, {
        wc_lang.Parameter: parameters,
        wc_lang.Species: all_species,
        wc_lang.Observable: {i.id: i for i in modifiers} if modifiers else {},
        wc_lang.Function: all_volumes,
    })
    assert error is None, str(error)
    print('---------------------------------------------------')
    return rate_law_expression, list(parameters.values())
