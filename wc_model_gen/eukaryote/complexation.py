""" Generator for macromolecular complexation submodel for eukaryotes
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-08-02
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import Bio.Alphabet
import Bio.Seq
import collections
import conv_opt
import wc_model_gen.utils as utils
import scipy.constants
import wc_kb
import wc_lang
import wc_model_gen


class ComplexationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for macromolecular complexation submodel

        Options:
        * amino_acid_id_conversion (:obj:`dict`): a dictionary with amino acid standard ids
            as keys and amino acid metabolite ids as values
        * codon_table (:obj:`dict`): a dictionary with protein subunit id as key and 
            NCBI identifier for translation table as value, the default is 1 (standard table) 
            for all protein
        * cds (:obj:`bool`): True indicates the sequences of protein subunits are complete CDS,
            the default if True    
        * beta (:obj:`float`, optional): ratio of Michaelis-Menten constant 
            to substrate concentration (Km/[S]) for use when estimating 
            Km values, the default value is 1
        * estimate_steady_state (:obj:`bool`): if True, the initial concentrations of complexes 
            and free-pool subunits will be estimated by assuming they are at steady state,
            the default if True          
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        if 'amino_acid_id_conversion' not in options:
            raise ValueError('The dictionary amino_acid_id_conversion has not been provided')
        else:    
            amino_acid_id_conversion = options['amino_acid_id_conversion']

        codon_table = options.get('codon_table', 1)
        options['codon_table'] = codon_table

        cds = options.get('cds', True)
        options['cds'] = cds    

        beta = options.get('beta', 1.)
        options['beta'] = beta

        estimate_steady_state = options.get('estimate_steady_state', True)
        options['estimate_steady_state'] = estimate_steady_state

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell

        amino_acid_id_conversion = self.options['amino_acid_id_conversion']
        codon_table = self.options['codon_table']
        cds = self.options['cds']
        
        print('Start generating complexation submodel...')
        assembly_rxn_no = 0
        disassembly_rxn_no = 0
        self._subunit_participation = collections.defaultdict(dict)
        for compl in cell.species_types.get(__type=wc_kb.core.ComplexSpeciesType):
            
            model_compl_species_type = model.species_types.get_one(id=compl.id)            
            
            for model_compl_species in model_compl_species_type.species:

                compl_compartment = model_compl_species.compartment
                
                if all(model.species_types.get_one(id=subunit.species_type.id).species.get_one(
                    compartment=compl_compartment)!=None for subunit in compl.subunits):
                    
                    # Generate complexation reactions
                    model_rxn = model.reactions.create(
                        submodel=self.submodel,
                        id='{}_association_in_{}'.format(compl.id, compl_compartment.id),
                        name='Complexation of {} in {}'.format(compl.id, compl_compartment.name),
                        reversible=False)

                    for subunit in compl.subunits:                        
                        model_subunit_species = model.species_types.get_one(
                            id=subunit.species_type.id).species.get_one(compartment=compl_compartment)
                        subunit_coefficient = subunit.coefficient if subunit.coefficient else 1
                        model_rxn.participants.add(
                            model_subunit_species.species_coefficients.get_or_create(
                            coefficient=-subunit_coefficient))
                        self._subunit_participation[model_subunit_species][model_compl_species] = subunit_coefficient

                    model_rxn.participants.add(
                        model_compl_species.species_coefficients.get_or_create(coefficient=1))

                    assembly_rxn_no += 1

                    # Generate reactions that lump complex dissociation and subunit degradation
                    for compl_subunit in compl.subunits:

                        if type(compl_subunit.species_type) == wc_kb.eukaryote.ProteinSpeciesType:                            
                            
                            model_rxn = model.reactions.create(
                                submodel=self.submodel,
                                id='{}_dissociation_in_{}_degradation_{}'.format(
                                    compl.id, compl_compartment.id, compl_subunit.species_type.id),
                                name='Dissociation of {} in {} and degradation of {}'.format(
                                    compl.id, compl_compartment.name, compl_subunit.species_type.id),
                                reversible=False)

                            model_rxn.participants.add(
                                model_compl_species.species_coefficients.get_or_create(
                                coefficient=-1))

                            disassembly_rxn_no += 1

                            for subunit in compl.subunits:
                                
                                model_subunit_species = model.species_types.get_one(
                                    id=subunit.species_type.id).species.get_one(compartment=compl_compartment)
                                subunit_coefficient = subunit.coefficient if subunit.coefficient else 1
                                
                                if subunit.species_type==compl_subunit.species_type:

                                    if codon_table == 1:
                                        codon_id = 1
                                    else:
                                        codon_id = codon_table[subunit.species_type.id]    
                                    
                                    protein_seq = ''.join(i for i in subunit.species_type.get_seq(table=codon_id, cds=cds) if i!='*')                                    
                                    aa_content = {}
                                    for aa in protein_seq:
                                        aa_id = amino_acid_id_conversion[aa]
                                        if aa_id not in aa_content:
                                            aa_content[aa_id] = 1
                                        else:
                                            aa_content[aa_id] += 1

                                    if compl_compartment.id == 'm':
                                        degradation_comp = model.compartments.get_one(id='m')
                                    else:
                                        degradation_comp = model.compartments.get_one(id='l')
                                    for aa_id, aa_count in aa_content.items():
                                        model_aa = model.species_types.get_one(id=aa_id).species.get_or_create(
                                            model=model, compartment=degradation_comp)
                                        model_rxn.participants.add(
                                            model_aa.species_coefficients.get_or_create(
                                            coefficient=aa_count))        

                                    h2o = model.species_types.get_one(id='h2o').species.get_one(
                                        compartment=degradation_comp)
                                    model_rxn.participants.add(
                                        h2o.species_coefficients.get_or_create(
                                        coefficient=-(sum(aa_content.values())-1)))                                   

                                    if subunit_coefficient > 1:
                                        model_rxn.participants.add(
                                            model_subunit_species.species_coefficients.get_or_create(
                                            coefficient=subunit_coefficient-1))

                                else:
                                    model_rxn.participants.add(
                                        model_subunit_species.species_coefficients.get_or_create(
                                        coefficient=subunit_coefficient))                              
        print('{} reactions of complex assembly and {} reactions of complex dissociation have been generated'.format(
            assembly_rxn_no, disassembly_rxn_no))

    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """
        model = self.model
        rate_law_no = 0
        for reaction in self.submodel.reactions:

            if 'association' in reaction.id:
                rate_law_exp, parameters = utils.gen_michaelis_menten_like_rate_law(
                    model, reaction)               

            else:
                diss_k_cat = model.parameters.create(id='k_cat_{}'.format(reaction.id),
                    type=wc_ontology['WC:k_cat'],
                    units=unit_registry.parse_units('s^-1'))
                
                complex_st_id = reaction.id[:reaction.id.index('_dissociation')]
                compl_compartment_id = reaction.id[reaction.id.index('_in_') + 4 : 
                    reaction.id.index('_degradation')]
                complex_species = model.species_types.get_one(id=complex_st_id).species.get_one(
                    compartment=model.compartments.get_one(id=compl_compartment_id))
                
                expression = '{} * {}'.format(diss_k_cat.id, complex_species.id)

                rate_law_exp, error = wc_lang.RateLawExpression.deserialize(expression, {
                    wc_lang.Parameter: {diss_k_cat.id: diss_k_cat},
                    wc_lang.Species: {complex_species.id: complex_species},
                })
                assert error is None, str(error)

            rate_law = model.rate_laws.create(
                    direction=wc_lang.RateLawDirection.forward,
                    type=None,
                    expression=rate_law_exp,
                    reaction=reaction,
                    )
            rate_law.id = rate_law.gen_id() 
            rate_law_no += 1
        
        print('{} rate laws for complex assembly and dissociation have been generated'.format(
            rate_law_no))

    def calibrate_submodel(self):
        """ Calibrate the submodel using data in the KB """        
        model = self.model        
        cell = self.knowledge_base.cell
        
        beta = self.options['beta']

        Avogadro = model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))

        # Calibrate dissociation constants
        for reaction in self.submodel.reactions:
            
            if 'dissociation' in reaction.id:
                complex_st_id = reaction.id[:reaction.id.index('_dissociation')]
                compl_compartment_id = reaction.id[reaction.id.index('_in_') + 4 : 
                    reaction.id.index('_degradation')]
                compl_compartment = model.compartments.get_one(id=compl_compartment_id)
                
                degraded_subunit_st_id = reaction.id[reaction.id.index('degradation_') + 12:]
                    
                degraded_subunit_hlife = cell.species_types.get_one(
                    id=degraded_subunit_st_id).properties.get_one(
                    property='half-life').get_value()
                degraded_subunit_species = model.species_types.get_one(
                    id=degraded_subunit_st_id).species.get_one(compartment=compl_compartment)
                degraded_subunit_stoic = model.reactions.get_one(id='{}_association_in_{}'.format(
                    complex_st_id, compl_compartment_id)).participants.get_one(
                    species=degraded_subunit_species).coefficient                         
                diss_k_cat = model.parameters.get_one(id='k_cat_{}'.format(reaction.id))
                diss_k_cat.value = - degraded_subunit_stoic / degraded_subunit_hlife

        # Calibrate the steady-state ('initial') concentrations of complex species by solving a system of linear equations
        if self.options['estimate_steady_state']:
            self.determine_steady_state_concentration(self._subunit_participation)    
                
        # Calibrate the parameter values for the rate laws of complex association reactions 
        for reaction in self.submodel.reactions:            
                        
            if 'association' in reaction.id:
                complex_st_id = reaction.id[:reaction.id.index('_association')]
                compl_compartment_id = reaction.id[reaction.id.index('_in_') + 4:]
                compl_compartment = model.compartments.get_one(id=compl_compartment_id)

                for param in reaction.rate_laws[0].expression.parameters:
                    if 'K_m_' in param.id:
                        species = model.species_types.get_one(
                            id=param.id.split('_')[-1]).species.get_one(
                            compartment=compl_compartment)
                        if species.distribution_init_concentration:
                            if species.distribution_init_concentration.mean:    
                                param.value = beta * species.distribution_init_concentration.mean \
                                    / Avogadro.value / species.compartment.init_volume.mean
                                param.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                                    beta, species.species_type.name, compl_compartment.name)    
                            else:
                                param.value = 1e-05
                                param.comments = 'The value was assigned to 1e-05 because the concentration of {} in {} was zero'.format(
                                    species.species_type.name, compl_compartment.name)        
                    elif 'k_cat_' in param.id:
                        param.value = 2e06
                        param.comments = 'The rate constant for bimolecular protein-protein association was used '\
                            'so that the simulated rate of complex assembly will be within the higher range'
                        param.references.append(wc_lang.Reference(
                            model=model,
                            title='Kinetics of protein-protein association explained by Brownian dynamics computer simulation',
                            author='Scott H Northrup, HP Erickson',        
                            year=1992,
                            type=wc_ontology['WC:article'],
                            publication='PNAS',        
                            volume='89',
                            issue='8',
                            pages='3338-3342'))        

        print('Complexation submodel has been generated')

    def determine_steady_state_concentration(self, subunit_participation):
        """ Use linear optimization to estimate the initial concentrations of complex species
            by assuming the system is at a state where the total amount of complex spexies is maximal. 
            The initial concentration of protein subunits will also be updated accordingly.

            Args:
                subunit_participation (:obj:`dict`): A nested dictionary with protein species as
                    keys, and dictionaries with key/value pairs of participated complex species and 
                    subunit coefficient as values
        """

        model = self.model

        all_variables = list(set([i.id for i in subunit_participation.keys()] + 
            [i.id for k,v in subunit_participation.items() for i in v.keys()]))
        
        opt_model = conv_opt.Model()
        
        variable_objs = {}
        for variable in all_variables:
            variable_objs[variable] = conv_opt.Variable(
                name=variable, type=conv_opt.VariableType.continuous, lower_bound=0)
            opt_model.variables.append(variable_objs[variable])
            if model.species.get_one(id=variable).species_type.type == wc_ontology['WC:pseudo_species']: 
                opt_model.objective_terms.append(
                    conv_opt.LinearTerm(variable_objs[variable], 1.))       
        opt_model.objective_direction = conv_opt.ObjectiveDirection.maximize      
        
        for subunit, participation in subunit_participation.items():            
            # Populate equation for subunit mass conservation
            constraint_coefs = [conv_opt.LinearTerm(variable_objs[subunit.id], 1)]                  
            for compl, coeff in participation.items():
                constraint_coefs.append(conv_opt.LinearTerm(variable_objs[compl.id], coeff))            
            if subunit.distribution_init_concentration: 
                total_mass = subunit.distribution_init_concentration.mean
            else:
                total_mass = 0.    
            constraint = conv_opt.Constraint(constraint_coefs, upper_bound=total_mass, lower_bound=total_mass)
            opt_model.constraints.append(constraint)
            
            # Populate inequation for steady-state of subunit concentration
            total_assoc_rate = 0
            constraint_coefs = []
            for compl, coeff in participation.items():
                
                total_diss_constant = 0
                
                dissociation_reactions = [i for i in model.reactions if '{}_dissociation_in_{}'.format(
                    compl.species_type.id, compl.compartment.id) in i.id]
                
                for reaction in dissociation_reactions:
                    if any(i.species==subunit for i in reaction.participants):
                        total_diss_constant += reaction.participants.get_one(species=subunit).coefficient * \
                            model.parameters.get_one(id='k_cat_{}'.format(reaction.id)).value
                
                constraint_coefs.append(conv_opt.LinearTerm(variable_objs[compl.id], total_diss_constant))   
                
                association_reaction = model.reactions.get_one(id='{}_association_in_{}'.format(
                    compl.species_type.id, compl.compartment.id))
                if all(i.distribution_init_concentration.mean > 0. for i in association_reaction.get_reactants() if i.distribution_init_concentration):
                    total_assoc_rate += 0.5 ** len([i for i in association_reaction.rate_laws[0].expression.parameters if 'K_m_' in i.id])    

            constraint = conv_opt.Constraint(constraint_coefs, upper_bound=2e06 * total_assoc_rate, lower_bound=0)
            opt_model.constraints.append(constraint)       
            
        # Solve the model and assign the results to species concentrations
        options = conv_opt.SolveOptions(solver=conv_opt.Solver.cplex, 
            presolve=conv_opt.Presolve.off)
        result = opt_model.solve(options)
        
        if result.status_code != conv_opt.StatusCode.optimal:
            raise Exception(result.status_message)
        
        primals = result.primals
        
        for ind, sol in enumerate(primals):               
            model_species = model.species.get_one(id=all_variables[ind])
            species_init_conc = model.distribution_init_concentrations.get_one(species=model_species)

            if species_init_conc:
                species_init_conc.mean = sol if sol > 0. else 0.
                species_init_conc.comments += '; Initial value was adjusted assuming the free pool ' + \
                        'is at steady state with its amount in macromolecular complexes'            

            else:
                conc_model = model.distribution_init_concentrations.create(
                    species=model_species,
                    mean=sol if sol > 0. else 0.,
                    units=unit_registry.parse_units('molecule'),
                    comments='Initial value was determined assuming the free pool ' + \
                        'is at steady state with its amount in macromolecular complexes'
                    )
                conc_model.id = conc_model.gen_id()
