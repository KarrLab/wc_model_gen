""" Generator for protein degradation submodel for eukaryotes
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-06-11
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.units import unit_registry
import math
import numpy
import scipy.constants
import wc_kb
import wc_lang
import wc_model_gen
import wc_model_gen.global_vars as gvar
import wc_model_gen.utils as utils


class ProteinDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for protein degradation submodel

        Options:
        * compartment_proteasomes (:obj:`dict`): a dictionary with compartment id
            as the key and a list of the names of proteasome complexes that degrade the protein
            species in the compartments as value
        * amino_acid_id_conversion (:obj:`dict`): a dictionary with amino acid standard ids
            as keys and amino acid metabolite ids as values
        * codon_table (:obj:`dict`): a dictionary with protein id as key and 
            NCBI identifier for translation table as value, the default is 1 (standard table) 
            for all protein
        * cds (:obj:`bool`): True indicates the sequences of protein are complete CDS,
            the default is True     
        * beta (:obj:`float`, optional): ratio of Michaelis-Menten constant 
            to substrate concentration (Km/[S]) for use when estimating 
            Km values, the default value is 1      
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        if 'compartment_proteasomes' not in options:
            raise ValueError('The dictionary of compartment:proteasomes has not been provided')
        else:    
            compartment_proteasomes = options['compartment_proteasomes']

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

        compartment_proteasomes = self.options['compartment_proteasomes']
        amino_acid_id_conversion = self.options['amino_acid_id_conversion']
        codon_table = self.options['codon_table']
        cds = self.options['cds']
                
        print('Start generating protein degradation submodel...')
        protein_kbs = cell.species_types.get(__type=wc_kb.eukaryote.ProteinSpeciesType)
        rxn_no = 0
        self._rxn_species_modifier = {}
        for protein_kb in protein_kbs:

            aa_content = {}

            if protein_kb.id in gvar.protein_aa_usage:                                        
                for aa, aa_id in amino_acid_id_conversion.items():
                    if gvar.protein_aa_usage[protein_kb.id][aa]:
                        aa_content[aa_id] = gvar.protein_aa_usage[protein_kb.id][aa]
            else:
                gvar.protein_aa_usage[protein_kb.id] = {i:0 for i in list(amino_acid_id_conversion.keys())}
                if codon_table == 1:
                    codon_id = 1
                else:
                    codon_id = codon_table[protein_kb.id]                                        
                raw_seq = protein_kb.get_seq(table=codon_id, cds=cds)                                            
                protein_seq = ''.join(i for i in raw_seq if i!='*')
                for aa in protein_seq:
                    aa_id = amino_acid_id_conversion[aa]
                    if aa_id not in aa_content:
                        aa_content[aa_id] = 1
                        gvar.protein_aa_usage[protein_kb.id][aa] = 1
                    else:
                        aa_content[aa_id] += 1
                        gvar.protein_aa_usage[protein_kb.id][aa] += 1
                gvar.protein_aa_usage[protein_kb.id]['*'] = raw_seq.count('*')
                gvar.protein_aa_usage[protein_kb.id]['len'] = len(protein_seq)
            
            protein_model = model.species_types.get_one(id=protein_kb.id)
            
            for protein_sp in protein_model.species:

                model_rxn = model.reactions.create(
                    submodel=self.submodel,
                    id='{}_{}_degradation'.format(protein_kb.id, protein_sp.compartment.id),
                    name='Degradation of {} of {}'.format(protein_kb.id, protein_sp.compartment.name),
                    reversible=False)

                self._rxn_species_modifier[model_rxn.id] = (protein_sp,)

                model_rxn.participants.add(
                    protein_sp.species_coefficients.get_or_create(coefficient=-1))

                if protein_sp.compartment.id not in compartment_proteasomes:
                    degradation_comp = model.compartments.get_one(id='l')
                else:
                    degradation_comp = protein_sp.compartment
                for aa_id, aa_count in aa_content.items():
                    model_aa = model.species_types.get_one(id=aa_id).species.get_or_create(
                        model=model, compartment=degradation_comp)
                    model_aa.id = model_aa.gen_id()
                    model_rxn.participants.add(
                        model_aa.species_coefficients.get_or_create(
                        coefficient=aa_count))        

                h2o = model.species_types.get_one(id='h2o').species.get_or_create(
                    compartment=degradation_comp)
                h2o.id = h2o.gen_id()
                model_rxn.participants.add(
                    h2o.species_coefficients.get_or_create(
                    coefficient=-(sum(aa_content.values())-1)))

                rxn_no += 1
        
        print('{} protein degradation reactions have been generated'.format(rxn_no))
            
    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """
        model = self.model
        compartment_proteasomes = self.options['compartment_proteasomes']
        
        rate_law_no = 0        
        for reaction in self.submodel.reactions:
            
            protein_compartment_id = self._rxn_species_modifier[reaction.id][0].compartment.id
            if protein_compartment_id not in compartment_proteasomes:
                degradation_comp = model.compartments.get_one(id='l')
            else:
                degradation_comp = model.compartments.get_one(id=protein_compartment_id)

            proteasomes = compartment_proteasomes[degradation_comp.id]
            if len(proteasomes) == 1:
                modifier = model.species_types.get_one(name=proteasomes[0]).species.get_one(
                    compartment=degradation_comp)
            else:
                if model.observables.get_one(id='total_proteasomes_{}'.format(degradation_comp.id)):
                    modifier = model.observables.get_one(id='total_proteasomes_{}'.format(
                        degradation_comp.id))
                else:
                    modifier_species = {}
                    
                    for proteasome in proteasomes:
                        proteasome_species = model.species_types.get_one(name=proteasome).species.get_one(
                            compartment=degradation_comp)
                        modifier_species[proteasome_species.id] = proteasome_species
                    
                    proteasome_total_exp, error = wc_lang.ObservableExpression.deserialize(
                        ' + '.join(modifier_species.keys()),
                        {wc_lang.Species: modifier_species})            
                    assert error is None, str(error)                
                    
                    modifier = model.observables.create(
                        id='total_proteasomes_{}'.format(degradation_comp.id), 
                        name='total proteasomes in {}'.format(degradation_comp.name), 
                        units=unit_registry.parse_units('molecule'), 
                        expression=proteasome_total_exp)

            self._rxn_species_modifier[reaction.id] += (modifier,)            

            h2o_species = model.species_types.get_one(id='h2o').species.get_one(
                compartment=degradation_comp)

            rate_law_exp, parameters = utils.gen_michaelis_menten_like_rate_law(
                self.model, reaction, modifiers=[modifier], exclude_substrates=[h2o_species])
            self.model.parameters += parameters

            rate_law = self.model.rate_laws.create(
                direction=wc_lang.RateLawDirection.forward,
                type=None,
                expression=rate_law_exp,
                reaction=reaction,
                )
            rate_law.id = rate_law.gen_id()
            rate_law_no += 1

        print('{} rate laws for protein degradation have been generated'.format(rate_law_no))     

    def calibrate_submodel(self):
        """ Calibrate the submodel using data in the KB """
        
        model = self.model        
        cell = self.knowledge_base.cell
        compartment_proteasomes = self.options['compartment_proteasomes']
        
        beta = self.options.get('beta')

        Avogadro = self.model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))       

        undetermined_model_kcat = []
        undetermined_model_km = []
        determined_kcat = []
        determined_km = []
        for reaction in self.submodel.reactions:

            init_species_counts = {}
            compartment_volumes = {}            
        
            modifier = self._rxn_species_modifier[reaction.id][1]
            modifier_total = 0      
            if type(modifier) == wc_lang.Species:
                init_species_counts[modifier.gen_id()] = modifier.distribution_init_concentration.mean
                modifier_total += modifier.distribution_init_concentration.mean
            else:
                for species in modifier.expression.species:
                    init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean
                    modifier_total += species.distribution_init_concentration.mean    
        
            protein_sp = self._rxn_species_modifier[reaction.id][0]
            mean_concentration = protein_sp.distribution_init_concentration.mean           
            half_life = cell.species_types.get_one(
                id=protein_sp.species_type.id).properties.get_one(property='half-life').get_value()
            
            average_rate = utils.calc_avg_deg_rate(mean_concentration, half_life)

            compartment_volumes[protein_sp.compartment.id] = protein_sp.compartment.init_volume.mean * \
                    protein_sp.compartment.init_density.value
            if protein_sp.compartment.id not in compartment_proteasomes:
                degradation_comp = model.compartments.get_one(id='l')                
            else:
                degradation_comp = model.compartments.get_one(id=protein_sp.compartment.id)
            
            model_Km = model.parameters.get_one(
                id='K_m_{}_{}'.format(reaction.id, protein_sp.species_type.id))
            if mean_concentration:
                model_Km.value = beta * mean_concentration \
                    / Avogadro.value / protein_sp.compartment.init_volume.mean
                if degradation_comp == protein_sp.compartment:    
                    model_Km.comments = 'The value was assumed to be ' +\
                        '{} times the concentration of {} in {}'.format(
                        beta, protein_sp.species_type.id, protein_sp.compartment.name)                     
                else:
                    model_Km.comments = 'The value was assumed to be ' +\
                        '{} times the concentration of {} in {} before transport to {}'.format(
                        beta, protein_sp.species_type.id, protein_sp.compartment.name, degradation_comp.name)
                determined_km.append(model_Km.value)
            else:
                undetermined_model_km.append(model_Km)            

            init_species_counts[protein_sp.gen_id()] = mean_concentration

            model_kcat = model.parameters.get_one(id='k_cat_{}'.format(reaction.id))

            if average_rate and modifier_total:
                model_kcat.value = 1.
                model_kcat.value = average_rate / reaction.rate_laws[0].expression._parsed_expression.eval({
                    wc_lang.Species: init_species_counts,
                    wc_lang.Compartment: compartment_volumes,
                })       
                determined_kcat.append(model_kcat.value)
            else:          
                undetermined_model_kcat.append(model_kcat)
        
        median_km = numpy.median(determined_km)
        for model_Km in undetermined_model_km:
            model_Km.value = median_km
            model_Km.comments = 'Set to the median value because protein concentration was zero'
        
        median_kcat = numpy.median(determined_kcat)
        for model_kcat in undetermined_model_kcat:
            model_kcat.value = median_kcat
            model_kcat.comments = 'Set to the median value because it could not be determined from data'       

        print('Protein degradation submodel has been generated')
