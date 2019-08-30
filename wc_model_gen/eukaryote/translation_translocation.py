""" Generator for translation and translocation submodel for eukaryotes
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-06-14
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.units import unit_registry
import wc_model_gen.utils as utils
import math
import numpy
import scipy.constants
import wc_kb
import wc_lang
import wc_model_gen
import wc_ontology


class TranslationTranslocationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for translation and translocation submodel 
        Options:
        * modifier_nucleus (:obj:`list`): a list of the names of complexes (ribosomes,
            initiation factors, elongation factors and release factors) involved in 
            translation in the nucleus
        * modifier_mitochondria (:obj:`list`): a list of the names of complexes (ribosomes,
            initiation factors, elongation factors and release factors) involved in 
            translation in the mitochondria
        * amino_acid_id_conversion (:obj:`dict`): a dictionary with amino acid standard ids
            as keys and amino acid metabolite ids as values         
        * beta (:obj:`float`, optional): ratio of Michaelis-Menten constant to substrate 
            concentration (Km/[S]) for use when estimating Km values, the default value is 1      
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        beta = options.get('beta', 1.)
        options['beta'] = beta

        if 'modifier_nucleus' not in options:
            raise ValueError('The list modifier_nucleus has not been provided')
        else:    
            modifier_nucleus = options['modifier_nucleus']

        if 'modifier_mitochondria' not in options:
            raise ValueError('The list modifier_mitochondria has not been provided')
        else:    
            modifier_mitochondria = options['modifier_mitochondria']

        if 'amino_acid_id_conversion' not in options:
            raise ValueError('The dictionary amino_acid_id_conversion has not been provided')
        else:    
            amino_acid_id_conversion = options['amino_acid_id_conversion']        

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        amino_acid_id_conversion = self.options.get('amino_acid_id_conversion')
        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')

        # Get tRNAs
        trna_kb = [i for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType) \
            if i.type==wc_kb.eukaryote.TranscriptType.tRna]
        trna_model = []
        for i in trna_kb:
            trna_model.append(model.species_types.get_one(id=i.id))        

        # Get species involved in reaction
        amino_acid_participants = list(amino_acid_id_conversion.values()) 
        other_metabolite_participants = ['gtp','gdp','pi']
        metabolites = {}
        for met in amino_acid_participants + other_metabolite_participants:
            met_species_type = model.species_types.get_one(id=met)
            metabolites[met] = {
                'n': met_species_type.species.get_one(compartment=nucleus),
                'm': met_species_type.species.get_one(compartment=mitochondrion)
                }            
        
        # Create reaction for each mRNA
        mrna_kbs = [i for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType) \
            if i.type==wc_kb.eukaryote.TranscriptType.mRna]
        
        for mrna_kb in mrna_kbs:
            
            mrna_kb_compartment_id = mrna_kb.species[0].compartment.id
            if mrna_kb_compartment_id == 'n':
                mrna_compartment = nucleus
            else:
                mrna_compartment = mitochondrion    
            
            mrna_model = model.species_types.get_one(id=mrna_kb.id).species.get_one(compartment=mrna_compartment)
            reaction = model.reactions.get_or_create(submodel=self.submodel, id='translation' + mrna_kb.id)
            reaction.name = 'translation of ' + mrna_kb.name
            mrna_seq = mrna_kb.get_seq()

            if mrna_compartment == mitochondrion:
                codon_table = 2
            else:
                codon_table = 1 
            protein_seq = mrna_kb.protein.get_seq(table=codon_table, cds=False).strip('*')

            # Adding participants to LHS
            reaction.participants.append(metabolites['gtp'][
                mrna_compartment.id].species_coefficients.get_or_create(coefficient=-(len(protein_seq)+2)))
            for aa_std, aa_met in amino_acid_id_conversion.items():
                if protein_seq.count(aa_std):
                    reaction.participants.append(metabolites[aa_met][
                        mrna_compartment.id].species_coefficients.get_or_create(coefficient=-protein_seq.count(aa_std)))            
            
            # Adding participants to RHS
            # TODO get and add protein species in all localised area
            reaction.participants.append(metabolites['gdp'][
                mrna_compartment.id].species_coefficients.get_or_create(coefficient=len(protein_seq)+2))
            reaction.participants.append(metabolites['pi'][
                mrna_compartment.id].species_coefficients.get_or_create(coefficient=len(protein_seq)*2))

            # Assign modifier
            polr_ids = rna_pol_pair[rna_kb.id]
            modifier_obs = []
            for polr_id in polr_ids: 
                complexes = model.species_types.get(name=polr_id)
                
                if not complexes:
                    raise ValueError('{} that catalyzes the transcription of {} cannot be found'.format(polr_id, rna_kb.id))
            
                else:                
                    observable = model.observables.get_one(
                        name='{} observable in {}'.format(polr_id, rna_compartment.name))

                if not observable:
                    
                    all_species = {}                             
                    
                    for compl_variant in complexes:
                        polr_species = compl_variant.species.get_one(compartment=rna_compartment)
                        if not polr_species:
                            raise ValueError('{} cannot be found in the {}'.format(polr_species, rna_compartment.name))
                        all_species[polr_species.gen_id()] = polr_species
                       
                    observable_expression, error = wc_lang.ObservableExpression.deserialize(
                        ' + '.join(list(all_species.keys())), {
                        wc_lang.Species: all_species,
                        })
                    assert error is None, str(error)

                    observable = model.observables.create(
                            name='{} observable in {}'.format(polr_id, rna_compartment.name),
                            expression=observable_expression)
                    observable.id = 'obs_{}'.format(len(model.observables))

                if observable not in modifier_obs:
                    modifier_obs.append(observable)

            if len(modifier_obs) == 1:
                self._transcription_modifier[reaction.name] = modifier_obs[0]
            else:
                all_obs = {obs.id: obs for obs in modifier_obs} 
                observable_expression, error = wc_lang.ObservableExpression.deserialize(
                    ' + '.join(list(all_obs.keys())), {
                    wc_lang.Observable: all_obs,
                    })
                assert error is None, str(error)

                observable = model.observables.create(
                            name='Combined RNA polymerase observable in {}'.format(rna_compartment.name),
                            expression=observable_expression)
                observable.id = 'obs_{}'.format(len(model.observables))

                self._transcription_modifier[reaction.name] = observable  

    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """                     

        modifier_nucleus = self.options.get('modifier_nucleus')
        modifier_mitochondria = self.options.get('modifier_mitochondria')
        for reaction in self.submodel.reactions:

            modifier = self._transcription_modifier[reaction.name]

            rate_law_exp, parameters = utils.gen_michaelis_menten_like_rate_law(
                self.model, reaction, modifiers=[modifier])
            self.model.parameters += parameters

            rate_law = self.model.rate_laws.create(
                direction=wc_lang.RateLawDirection.forward,
                type=None,
                expression=rate_law_exp,
                reaction=reaction,
                )
            rate_law.id = rate_law.gen_id()
                        
    def calibrate_submodel(self):
        """ Calibrate the submodel using data in the KB """
        
        model = self.model        
        cell = self.knowledge_base.cell
        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')

        beta = self.options.get('beta')

        Avogadro = self.model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))

        mean_doubling_time = model.parameters.get_one(id='mean_doubling_time').value       
        
        rnas_kb = cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType)
        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):

            init_species_counts = {}
        
            modifier = self._transcription_modifier[reaction.name]      
            for species in modifier.expression.species:
                init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean
            for observable in modifier.expression.observables:
                for species in observable.expression.species:
                    init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean    

            rna_kb_compartment_id = rna_kb.species[0].compartment.id
            if rna_kb_compartment_id == 'n':
                rna_compartment = nucleus
            else:
                rna_compartment = mitochondrion    
            
            rna_product = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=rna_compartment)           
            
            half_life = rna_kb.properties.get_one(property='half_life').get_value()
            mean_concentration = rna_product.distribution_init_concentration.mean         

            average_rate = utils.calc_avg_syn_rate(
                mean_concentration, half_life, mean_doubling_time)
            
            for species in reaction.get_reactants():
                
                init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean
                
                if model.parameters.get(id='K_m_{}_{}'.format(reaction.id, species.species_type.id)):
                    model_Km = model.parameters.get_one(
                        id='K_m_{}_{}'.format(reaction.id, species.species_type.id))
                    model_Km.value = beta * species.distribution_init_concentration.mean \
                        / Avogadro.value / species.compartment.init_volume.mean

            model_kcat = model.parameters.get_one(id='k_cat_{}'.format(reaction.id))
            model_kcat.value = 1.
            model_kcat.value = average_rate / reaction.rate_laws[0].expression._parsed_expression.eval({
                wc_lang.Species: init_species_counts,
                wc_lang.Compartment: {
                    rna_compartment.id: rna_compartment.init_volume.mean * rna_compartment.init_density.value},
                })