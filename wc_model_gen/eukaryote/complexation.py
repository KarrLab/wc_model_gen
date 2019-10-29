""" Generator for macromolecular complexation submodel for eukaryotes
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-08-02
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import wc_model_gen.global_vars as gvar
import Bio.Alphabet
import Bio.Seq
import collections
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
        * estimate_initial_state (:obj:`bool`): if True, the initial concentrations of complexes 
            and free-pool subunits will be estimated using linear programming,
            the default is True
        * greedy_step_size (:obj:`float`): the extent to which complex copy number is increased
            at each round of reaction selection during initial copy number estimation, 
            value should be higher than 0 and not more than 1.0, and the default value is 0.1    
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

        estimate_initial_state = options.get('estimate_initial_state', True)
        options['estimate_initial_state'] = estimate_initial_state

        greedy_step_size = options.get('greedy_step_size', 0.1)
        assert(greedy_step_size > 0 and greedy_step_size <= 1.)
        options['greedy_step_size'] = greedy_step_size

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
        self._maximum_possible_amount = {}
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
                    
                    self._maximum_possible_amount[model_compl_species.id] = []
                    for subunit in compl.subunits:                        
                        model_subunit_species = model.species_types.get_one(
                            id=subunit.species_type.id).species.get_one(compartment=compl_compartment)
                        subunit_coefficient = subunit.coefficient if subunit.coefficient else 1
                        model_rxn.participants.add(
                            model_subunit_species.species_coefficients.get_or_create(
                            coefficient=-subunit_coefficient))
                        if model_subunit_species.species_type.type != wc_ontology['WC:metabolite']:
                            self._maximum_possible_amount[model_compl_species.id].append(
                                model_subunit_species.distribution_init_concentration.mean / subunit_coefficient)

                    model_rxn.participants.add(
                        model_compl_species.species_coefficients.get_or_create(coefficient=1))

                    assembly_rxn_no += 1

                    # Generate reactions that lump complex dissociation and subunit degradation
                    for compl_subunit in compl.subunits:

                        if type(compl_subunit.species_type) == wc_kb.eukaryote.ProteinSpeciesType:                            
                            
                            model_rxn = model.reactions.create(
                                submodel=self.submodel,
                                id='{}_dissociation_in_{}_degrade_{}'.format(
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
                                    
                                    aa_content = {}

                                    if subunit.species_type.id in gvar.protein_aa_usage:                                        
                                        for aa, aa_id in amino_acid_id_conversion.items():
                                            if gvar.protein_aa_usage[subunit.species_type.id][aa]:
                                                aa_content[aa_id] = gvar.protein_aa_usage[subunit.species_type.id][aa]
                                    else:
                                        if codon_table == 1:
                                            codon_id = 1
                                        else:
                                            codon_id = codon_table[subunit.species_type.id]                                        
                                        protein_seq = ''.join(i for i in subunit.species_type.get_seq(table=codon_id, cds=cds) if i!='*')                                   
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
        
        self._maximum_possible_amount = {k:min(v) if v else 0 for k, v in self._maximum_possible_amount.items()}

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
                rate_law_exp.expression += ' * 2 ** {}'.format(len(reaction.get_reactants()))               

            else:
                diss_k_cat = model.parameters.create(id='k_cat_{}'.format(reaction.id),
                    type=wc_ontology['WC:k_cat'],
                    units=unit_registry.parse_units('molecule^-1 s^-1'))
                
                complex_st_id = reaction.id[:reaction.id.index('_dissociation')]
                compl_compartment_id = reaction.id[reaction.id.index('_in_') + 4 : 
                    reaction.id.index('_degrade')]
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
        self._effective_dissociation_constant = collections.defaultdict(float)
        for reaction in self.submodel.reactions:
            
            if 'dissociation' in reaction.id:
                complex_st_id = reaction.id[:reaction.id.index('_dissociation')]
                compl_compartment_id = reaction.id[reaction.id.index('_in_') + 4 : 
                    reaction.id.index('_degrade')]
                compl_compartment = model.compartments.get_one(id=compl_compartment_id)
                
                degraded_subunit_st_id = reaction.id[reaction.id.index('degrade_') + 8:]
                    
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
                
                model_compl_species = model.species_types.get_one(
                    id=complex_st_id).species.get_one(compartment=compl_compartment)
                self._effective_dissociation_constant[model_compl_species.id] += diss_k_cat.value

        # Estimate the initial concentrations of complex species
        if self.options['estimate_initial_state']:
            self.determine_initial_concentration()    
                
        # Calibrate the parameter values for the rate laws of complex association reactions 
        ref_kcat = wc_lang.Reference(
            model=model,
            title='Kinetics of protein-protein association explained by Brownian dynamics computer simulation',
            author='Scott H Northrup, HP Erickson',        
            year=1992,
            type=wc_ontology['WC:article'],
            publication='PNAS',        
            volume='89',
            issue='8',
            pages='3338-3342')
        ref_kcat.id = 'ref_'+str(len(model.references))
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
                                    beta, species.species_type.id, compl_compartment.name)    
                            else:
                                param.value = 1e-05
                                param.comments = 'The value was assigned to 1e-05 because the concentration of {} in {} was zero'.format(
                                    species.species_type.id, compl_compartment.name)        
                    elif 'k_cat_' in param.id:
                        param.value = 2e06
                        param.comments = 'The rate constant for bimolecular protein-protein association was used '\
                            'so that the simulated rate of complex assembly will be within the higher range'
                        param.references.append(ref_kcat)        

        print('Complexation submodel has been generated')

    def determine_initial_concentration(self):
        """ Estimate the initial concentrations of complex species using the following steps: 

            1) For each complex species, calculate the maximum possible copy number by taking the minimum
               of the availability of each subunit, which is determined as the ratio of
               subunit copy number to its stoichiometric coefficient in the complex.
            2) Arrange complexation reactions in decreasing order of the effective 
               dissociation rate into a list.
            3) For each reaction in the list, increase the copy number of complex 
               by either the minimum of all current subunit availability or the maximum possible
               copy number calculated in step 1 multiplied by the greedy_step_size, whichever is less. 
               The copy number of subunits are adjusted accordingly. If the copy
               number of any subunits reaches zero, the reaction is removed from the list.  
            4) Repeat step 3 until the list is empty.
        """

        model = self.model
        greedy_step_size = self.options['greedy_step_size']
        
        complexation_reactions = {}
        for reaction in self.submodel.reactions:
            if 'association' in reaction.id:
                complex_species = reaction.get_products()[0]
                complexation_reactions[reaction] = self._effective_dissociation_constant[complex_species.id] * \
                    self._maximum_possible_amount[complex_species.id]

        sorted_reactions = sorted(complexation_reactions.items(), key=lambda kv: kv[1], reverse=True)
        sorted_reactions = [(v1, v2) for v1, v2 in sorted_reactions if v2]
        sorted_reactions = collections.OrderedDict(sorted_reactions)

        while sorted_reactions:
            for_removal = []
            for reaction in sorted_reactions:
                complex_species = reaction.get_products()[0]
                
                current_availability = []
                for participant in reaction.participants:
                    if participant.coefficient < 0 and participant.species.species_type.type != wc_ontology['WC:metabolite']:
                        current_availability.append(- participant.species.distribution_init_concentration.mean / participant.coefficient)
                
                flux = min(greedy_step_size*self._maximum_possible_amount[complex_species.id], min(current_availability))
                for participant in reaction.participants:
                    if participant.species.species_type.type != wc_ontology['WC:metabolite']:
                        species_init_conc = model.distribution_init_concentrations.get_one(species=participant.species)

                        if species_init_conc:
                            species_init_conc.mean += flux * participant.coefficient
                            
                        else:
                            conc_model = model.distribution_init_concentrations.create(
                                species=participant.species,
                                mean=flux * participant.coefficient,
                                units=unit_registry.parse_units('molecule'),
                                )
                            conc_model.id = conc_model.gen_id()

                    if participant.species.species_type.type == wc_ontology['WC:protein']:
                        if model.distribution_init_concentrations.get_one(species=participant.species).mean == 0.:
                            if reaction not in for_removal:
                                for_removal.append(reaction)
            
            for reaction in for_removal:
                del sorted_reactions[reaction]                            
