""" Utility methods for generating submodels

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-23
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import collections
import conv_opt
import copy
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


def test_metabolite_production(submodel, reaction_bounds, pseudo_reactions=None, 
    test_products=None, test_reactants=None):
    """ Test that an FBA metabolism submodel can produce each product component and 
        recycle each reactant component in each reaction in the objective function
        individually. For each product (reactant), a sink (source) reaction is added 
        as objective function to be maximized. If the solution returns a zero 
        objective function, that indicates the submodel cannot produce (recycle) the
        product (reactant).

        Args:
            submodel (:obj:`wc_lang.Submodel`): FBA metabolism submodel
            reaction_bounds (:obj:`dict`): a dictionary with reaction IDs as keys
                and tuples of (lower bound, upper bound) as values
            pseudo_reactions (:obj:`list` of :obj:`str`, optional): list of IDs of 
                pseudo-reactions that will be excluded from the optimization formulation
            test_products (:obj:`list` of :obj:`str`, optional): list of IDs of species 
                to be tested for producibility in the submodel; if nothing is provided, 
                the species will be extracted from the products in objective reactions
            test_reactants (:obj:`list` of :obj:`str`, optional): list of IDs of species 
                to be tested for recyclability in the submodel; if nothing is provided, 
                the species will be extracted from the reactants in objective reactions
            
        Returns:
            :obj:`list`: list of metabolite species IDs that cannot be produced 
                in the metabolism submodel
            :obj:`list`: list of metabolite species IDs that cannot be recycled 
                in the metabolism submodel         
    """
    if test_products:
        products = test_products
    else:    
        products = []
        for reaction in submodel.dfba_obj.expression.reactions:
            for product in reaction.get_products():
                if product.id not in products:
                    products.append(product.id)
    
    if test_reactants:
        reactants = test_reactants
    else:    
        reactants = []
        for reaction in submodel.dfba_obj.expression.reactions:
            for reactant in reaction.get_reactants():
                if reactant.id not in reactants:
                    reactants.append(reactant.id)

    if pseudo_reactions:
        excluded_rxns = pseudo_reactions
    else:
        excluded_rxns = []                    

    conv_model = conv_opt.Model(name='model')
    conv_variables = {}
    conv_metabolite_matrices = collections.defaultdict(list)
    for reaction in submodel.reactions:
        if reaction.id not in excluded_rxns:
            conv_variables[reaction.id] = conv_opt.Variable(
                name=reaction.id, type=conv_opt.VariableType.continuous,
                lower_bound=reaction_bounds[reaction.id][0], 
                upper_bound=reaction_bounds[reaction.id][1])
            conv_model.variables.append(conv_variables[reaction.id])
            for part in reaction.participants:
                conv_metabolite_matrices[part.species.id].append(
                    conv_opt.LinearTerm(conv_variables[reaction.id], 
                        part.coefficient))  

    for met_id, expression in conv_metabolite_matrices.items():
        conv_model.constraints.append(conv_opt.Constraint(expression, name=met_id, 
            upper_bound=0.0, lower_bound=0.0))                    
          
    conv_model.objective_direction = conv_opt.ObjectiveDirection.maximize

    options = conv_opt.SolveOptions(
        solver=conv_opt.Solver.cplex,
        presolve=conv_opt.Presolve.on,
        solver_options={
            'cplex': {
                'parameters': {
                    'emphasis': {
                        'numerical': 1,
                    },
                    'read': {
                        'scale': 1,
                    },
                },
            },
        })

    unproducibles = []
    for product in products:
        temp_model = copy.deepcopy(conv_model)        
        obj_reaction = conv_opt.Variable(name=product + '_test_reaction', 
                                        type=conv_opt.VariableType.continuous,
                                        lower_bound=0.0)
        temp_model.variables.append(obj_reaction)        
        
        expression = [i for i in temp_model.constraints if i.name==product]
        if expression:
            expression[0].terms.append(conv_opt.LinearTerm(obj_reaction, 1.))
        else:
            temp_model.constraints.append(conv_opt.Constraint(
                [conv_opt.LinearTerm(obj_reaction, 1.)], name=product, 
                upper_bound=0.0, lower_bound=0.0))
        
        temp_model.objective_terms = [conv_opt.LinearTerm(obj_reaction, 1.),]
        result = temp_model.solve()
        if result.value == 0.0:
            unproducibles.append(product)            

    unrecyclables = []
    for reactant in reactants:
        temp_model = copy.deepcopy(conv_model)
        obj_reaction = conv_opt.Variable(name=reactant + '_test_reaction', 
                                        type=conv_opt.VariableType.continuous,
                                        lower_bound=0.0)
        temp_model.variables.append(obj_reaction)
        
        expression = [i for i in temp_model.constraints if i.name==reactant]
        if expression:
            expression[0].terms.append(conv_opt.LinearTerm(obj_reaction, -1.))
        else:
            temp_model.constraints.append(conv_opt.Constraint(
                [conv_opt.LinearTerm(obj_reaction, -1.)], name=reactant, 
                upper_bound=0.0, lower_bound=0.0))
        
        temp_model.objective_terms = [conv_opt.LinearTerm(obj_reaction, 1.),]
        result = temp_model.solve()
        if result.value == 0.0:
            unrecyclables.append(reactant) 

    return unproducibles, unrecyclables                 


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
        type=None,
        units=unit_registry.parse_units(''))
    parameters[f.id] = f

    volume = activator.compartment.init_density.function_expressions[0].function
    functions[volume.id] = volume
    
    F_act = '((1 + {} / ({} * {} * {}) * {}) / (1 + {} / ({} * {} * {})))'.format(
        activator.id, Ka.id, avogadro.id, volume.id, f.id, activator.id, Ka.id, avogadro.id, volume.id)   

    return F_act, species, parameters, functions


def gen_michaelis_menten_like_rate_law(model, reaction, modifiers=None, modifier_reactants=None, exclude_substrates=None):
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
            modifiers (:obj:`list` of :obj:`wc_lang.Observable` or :obj:`wc_lang.Species`): 
                list of observables (each of which evaluates to the total concentration of all enzymes 
                that catalyze the same intermediate step in the reaction) or enzyme species
                (that catalyze the reaction) 
            modifier_reactants (:obj:`list` of :obj:`wc_lang.Species`): list of species 
                in modifiers that should be included as reactants in the rate law 
            exclude_substrates (:obj:`list` of :obj:`wc_lang.Species`): list of reactant species 
                that would be excluded from the rate law       

        Returns:
                :obj:`wc_lang.RateLawExpression`: rate law
                :obj:`list` of :obj:`wc_lang.Parameter`: list of parameters in the rate law     
    """
    modifier_species = []
    all_species = {}
    all_volumes = {}
    all_observables = {}
    all_parameters = {}
    if modifiers:
        for modifier in modifiers:
            if type(modifier) == wc_lang.Observable:
                all_observables[modifier.id] = modifier
                for species in modifier.expression.species:
                    modifier_species.append(species)                    
            elif type(modifier) == wc_lang.Species:
                modifier_species.append(modifier)
                all_species[modifier.gen_id()] = modifier
            else:
                raise TypeError('The modifiers contain element(s) that is not an observable or a species')          

    if modifier_reactants is None:
        additional_reactants = []
    else:
        additional_reactants = modifier_reactants

    if exclude_substrates:
        excluded_reactants = exclude_substrates
    else:
        excluded_reactants = []    

    avogadro = model.parameters.get_or_create(
        id='Avogadro',
        type=None,
        value=scipy.constants.Avogadro,
        units=unit_registry.parse_units('molecule mol^-1'))
    all_parameters[avogadro.id] = avogadro

    model_k_cat = model.parameters.get_or_create(id='k_cat_{}'.format(reaction.id),
                                                 type=wc_ontology['WC:k_cat'],
                                                 units=unit_registry.parse_units('s^-1{}'.format(
                                                (' * molecule^{{-{}}}'.format(len(modifiers))) if modifiers else '')))
    all_parameters[model_k_cat.id] = model_k_cat

    expression_terms = []    
    for species in reaction.get_reactants():

        if (species not in modifier_species or species in additional_reactants) and species not in excluded_reactants:

            all_species[species.gen_id()] = species

            model_k_m = model.parameters.create(id='K_m_{}_{}'.format(reaction.id, species.species_type.id),
                                                type=wc_ontology['WC:K_m'],
                                                units=unit_registry.parse_units('M'))
            all_parameters[model_k_m.id] = model_k_m

            volume = species.compartment.init_density.function_expressions[0].function
            all_volumes[volume.id] = volume

            expression_terms.append('({} / ({} + {} * {} * {}))'.format(species.gen_id(),
                                                                        species.gen_id(),
                                                                        model_k_m.id, avogadro.id,
                                                                        volume.id))

    expression = '{}{}{}'.format(
        model_k_cat.id,
        (' * {}'.format(' * '.join([i.id if type(i)==wc_lang.Observable else i.gen_id() \
            for i in modifiers]))) if modifiers else '',
        (' * {}'.format(' * '.join(expression_terms))) if expression_terms else '')
    
    rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, {
        wc_lang.Parameter: all_parameters,
        wc_lang.Species: all_species,
        wc_lang.Observable: all_observables,
        wc_lang.Function: all_volumes,
    })
    assert error is None, str(error)

    return rate_law_expression, list(all_parameters.values())


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

def gen_response_functions(model, beta, reaction_id, reaction_class, compartment, reaction_factors):
        """ Generate a list of response function expression string for each factor or 
            group of factors (F) in the form of:
                       
                        F/(Km + F)

        Args:
            model (:obj:`wc_lang.Model`): model
            beta (:obj:`float`): ratio of Michaelis-Menten constant to substrate 
                concentration (Km/[S]) for use when estimating Km values
            reaction_id (:obj:`str`): identifier of reaction whose rate law will use the function expressions
            reaction_class (:obj:`str`): generic class that the reaction belongs to which 
                shares the same observables in their rate laws, e.g. 'translation_initiation' for each gene that 
                shares the same initiation factors 
            compartment (:obj:`wc_lang.Compartment`): compartment where the reaction occurs
            reaction_factors (:obj:`list` of `list`): list of lists of the ID or name of
                (initiation/elongation/translocation) factors, grouped based on similar functions or classes,
                e.g. [['factor1 variant1', 'factor1 variant2'], ['factor2']]
            
        Returns:
            :obj:`list`: list of strings of response function expression for each factor/group of factors
            :obj:`dict`: IDs of species (keys) and their species objects (values)
            :obj:`dict`: IDs of parameters (keys) and their parameter objects (values)
            :obj:`dict`: IDs of volume density functions (keys) and their function objects (values)
            :obj:`dict`: IDs of observables (keys) and their observable objects (values)
        """
        
        all_species = {}
        all_parameters = {}
        all_volumes = {}
        all_observables = {}

        Avogadro = model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))
        all_parameters[Avogadro.id] = Avogadro

        volume = compartment.init_density.function_expressions[0].function
        all_volumes[volume.id] = volume

        factor_exp = []
        for factors in reaction_factors:
            
            if len(factors) == 1:

                factor_species_type = model.species_types.get_one(name=factors[0]) 
                if not factor_species_type:
                    factor_species_type = model.species_types.get_one(id=factors[0])
                factor_species = factor_species_type.species.get_one(compartment=compartment)                
                all_species[factor_species.gen_id()] = factor_species

                model_k_m = model.parameters.create(
                    id='K_m_{}_{}'.format(reaction_id, factor_species.species_type.id),
                    value = beta * factor_species.distribution_init_concentration.mean \
                        / Avogadro.value / compartment.init_volume.mean,
                    type=wc_ontology['WC:K_m'],
                    units=unit_registry.parse_units('M'),
                    comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                        beta, factor_species_type.id, compartment.name)
                    )
                if model_k_m.value == 0.:
                    model_k_m.value = 1e-05
                    model_k_m.comments = 'The value was assigned to 1e-05 because the concentration of ' \
                        '{} in {} was zero'.format(factor_species_type.id, compartment.name)
                all_parameters[model_k_m.id] = model_k_m                    

                factor_exp.append('({} / ({} + {} * {} * {}))'.format(
                    factor_species.gen_id(),
                    factor_species.gen_id(),
                    model_k_m.id, 
                    Avogadro.id,
                    volume.id))

            else:
                
                obs_exp = []
                obs_total = 0                    
                for factor in factors:                  
                    factor_species_type = model.species_types.get_one(name=factor)
                    if not factor_species_type:
                        factor_species_type = model.species_types.get_one(id=factor)
                    factor_species = factor_species_type.species.get_one(compartment=compartment)
                    all_species[factor_species.gen_id()] = factor_species
                    obs_exp.append(factor_species.gen_id())
                    obs_total += factor_species.distribution_init_concentration.mean
                    obs_exp_string = ' + '.join(sorted(obs_exp))
                
                if not any(i.expression.expression==obs_exp_string for i in model.observables):
                    observable_exp, error = wc_lang.ObservableExpression.deserialize(
                    obs_exp_string,
                    {wc_lang.Species: all_species})            
                    assert error is None, str(error)                
                    
                    n = len(model.observables.get(name='factor for {} in {}'.format(
                        reaction_class, compartment.name)))

                    factor_observable = model.observables.get_or_create(
                        id='{}_factors_{}_{}'.format(reaction_class, compartment.id, n+1), 
                        name='factor for {} in {}'.format(reaction_class, compartment.name), 
                        units=unit_registry.parse_units('molecule'), 
                        expression=observable_exp)
                    all_observables[factor_observable.id] = factor_observable
                else:
                    factor_observable = [i for i in model.observables if i.expression.expression==obs_exp_string][0]
                    all_observables[factor_observable.id] = factor_observable

                model_k_m = model.parameters.create(
                    id='K_m_{}_{}'.format(reaction_class, factor_observable.id),
                    value = beta * obs_total / Avogadro.value / compartment.init_volume.mean,
                    type=wc_ontology['WC:K_m'],
                    units=unit_registry.parse_units('M'),
                    comments = 'The value was assumed to be {} times the value of {}'.format(
                        beta, factor_observable.id)  
                    )
                if model_k_m.value == 0.:
                    model_k_m.value = 1e-05
                    model_k_m.comments = 'The value was assigned to 1e-05 because the value of ' \
                        '{} was zero'.format(factor_observable.id)
                all_parameters[model_k_m.id] = model_k_m

                factor_exp.append('({} / ({} + {} * {} * {}))'.format(
                    factor_observable.id,
                    factor_observable.id,
                    model_k_m.id, 
                    Avogadro.id,
                    volume.id))

        return factor_exp, all_species, all_parameters, all_volumes, all_observables

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
