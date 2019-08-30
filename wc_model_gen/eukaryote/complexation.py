""" Generator for macromolecular complexation submodel for eukaryotes
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-08-02
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
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
        for compl in cell.species_types.get(__type=wc_kb.core.ComplexSpeciesType):
            
            model_compl_species_type = model.species_types.get_one(id=compl.id)            
            
            for model_compl_species in model_compl_species_type.species:

                compl_compartment = model_compl_species.compartment
                
                if all(model.species_types.get_one(id=subunit.species_type.id).species.get_one(
                    compartment=compl_compartment)!=None for subunit in compl.subunits):
                    
                    # Generate complexation reactions
                    model_rxn = model.reactions.create(
                        submodel=self.submodel,
                        id='complex_association_{}_{}'.format(compl.id, compl_compartment.id),
                        name='Complexation of {} in {}'.format(compl.id, compl_compartment.name),
                        reversible=False)

                    for subunit in compl.subunits:                        
                        model_subunit_species = model.species_types.get_one(
                            id=subunit.species_type.id).species.get_one(compartment=compl_compartment)
                        subunit_coefficient = subunit.coefficient if subunit.coefficient else 1
                        model_rxn.participants.add(
                            model_subunit_species.species_coefficients.get_or_create(
                            coefficient=-subunit_coefficient))

                    model_rxn.participants.add(
                        model_compl_species.species_coefficients.get_or_create(coefficient=1))

                    assembly_rxn_no += 1

                    # Generate reactions that lump complex dissociation and subunit degradation
                    for compl_subunit in compl.subunits:

                        if type(compl_subunit.species_type) == wc_kb.eukaryote.ProteinSpeciesType:                            
                            
                            model_rxn = model.reactions.create(
                                submodel=self.submodel,
                                id='{}_{}_dissociation_{}_degradation'.format(
                                    compl.id, compl_compartment.id, compl_subunit.species_type.id),
                                name='Dissociation of {} and degradation of {} in {}'.format(
                                    compl.id, compl_subunit.species_type.id, compl_compartment.name),
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
                                    
                                    protein_seq = subunit.species_type.get_seq(table=codon_id, cds=cds)
                                    aa_content = {}
                                    for aa in protein_seq:
                                        aa_id = amino_acid_id_conversion[aa]
                                        if aa_id not in aa_content:
                                            aa_content[aa_id] = 1
                                        else:
                                            aa_content[aa_id] += 1

                                    for aa_id, aa_count in aa_content.items():
                                        model_aa = model.species_types.get_one(id=aa_id).species.get_one(
                                            compartment=compl_compartment)
                                        model_rxn.participants.add(
                                            model_aa.species_coefficients.get_or_create(
                                            coefficient=aa_count))        

                                    h2o = model.species_types.get_one(id='h2o').species.get_one(
                                        compartment=compl_compartment)
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

            if 'complex_association_' in reaction.id:
                rate_law_exp, parameters = utils.gen_michaelis_menten_like_rate_law(
                    model, reaction)               

            else:
                diss_k_cat = model.parameters.create(id='k_cat_{}'.format(reaction.id),
                    type=wc_ontology['WC:k_cat'],
                    units=unit_registry.parse_units('s^-1'))
                
                complex_species = model.species_types.get_one(id=reaction.id.split('_')[0]).species.get_one(
                    compartment=model.compartments.get_one(id=reaction.id.split('_')[1]))
                
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
        
        print('{} rate laws for complex assembly and dissociation have been generated'.format(rate_law_no))

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

        for reaction in self.submodel.reactions:            
            
            if 'complex_association_' in reaction.id:

                compl_compartment = model.compartments.get_one(id=reaction.id.split('_')[-1])

                for param in reaction.rate_laws[0].expression.parameters:
                    if 'K_m_' in param.id:
                        species = model.species_types.get_one(
                            id=param.id.split('_')[-1]).species.get_one(
                            compartment=compl_compartment)
                        param.value = beta * species.distribution_init_concentration.mean \
                            / Avogadro.value / species.compartment.init_volume.mean
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
            
            else:

                compl_compartment = model.compartments.get_one(id=reaction.id.split('_')[1])

                degraded_subunit_hlife = cell.species_types.get_one(
                    id=reaction.id.split('_')[3]).properties.get_one(
                    property='half_life').get_value()

                degraded_subunit_species = model.species_types.get_one(
                    id=reaction.id.split('_')[3]).species.get_one(compartment=compl_compartment)
                degraded_subunit_stoic = model.reactions.get_one(id='complex_association_{}_{}'.format(
                    reaction.id.split('_')[0], compl_compartment.id)).participants.get_one(
                    species=degraded_subunit_species).coefficient                           

                diss_k_cat = model.parameters.get_one(id='k_cat_{}'.format(reaction.id))
                diss_k_cat.value = - degraded_subunit_stoic / degraded_subunit_hlife

        print('Complexation submodel has been generated')                                
