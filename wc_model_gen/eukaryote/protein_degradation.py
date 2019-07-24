""" Generator for protein degradation submodel for eukaryotes
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-06-11
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


class ProteinDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for protein degradation submodel

        Options:
        * protein_proteasomes (:obj:`dict`): a dictionary with protein species ID
            as the key and a list of the names of proteasome complexes that degrade the protein
            species as value 
        * beta (:obj:`float`, optional): ratio of Michaelis-Menten constant 
            to substrate concentration (Km/[S]) for use when estimating 
            Km values, the default value is 1      
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        beta = options.get('beta', 1.)
        options['beta'] = beta

        if 'protein_proteasomes' not in options:
            raise ValueError('The list protein_proteasomes has not been provided')
        else:    
            protein_proteasomes = options['protein_proteasomes']

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
                
        # Get species involved in reaction
        metabolic_participants = ['amp', 'cmp', 'gmp', 'ump', 'h2o', 'h']
        metabolites = {}
        for met in metabolic_participants:
            met_species_type = model.species_types.get_one(id=met)
            metabolites[met] = {
                'n': met_species_type.species.get_one(compartment=nucleus),
                'm': met_species_type.species.get_one(compartment=mitochondrion)
                }

        # Create reaction for each RNA and get exosome
        rna_exo_pair = self.options.get('rna_exo_pair')
        rna_kbs = cell.species_types.get(__type=wc_kb.eukaryote_schema.TranscriptSpeciesType)
        self._degradation_modifier = {}
        for rna_kb in rna_kbs:  

            rna_kb_compartment_id = rna_kb.species[0].compartment.id
            if rna_kb_compartment_id == 'n':
                rna_compartment = nucleus
            else:
                rna_compartment = mitochondrion    
            
            rna_model = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=rna_compartment)
            reaction = model.reactions.get_or_create(submodel=self.submodel, id='degradation_' + rna_kb.id)
            reaction.name = 'degradation of ' + rna_kb.name
            seq = rna_kb.get_seq()

            # Adding participants to LHS
            reaction.participants.append(rna_model.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.append(metabolites['h2o'][
                rna_compartment.id].species_coefficients.get_or_create(coefficient=-(len(seq)-1)))

            # Adding participants to RHS
            reaction.participants.append(metabolites['amp'][
                rna_compartment.id].species_coefficients.get_or_create(coefficient=seq.count('A')))
            reaction.participants.append(metabolites['cmp'][
                rna_compartment.id].species_coefficients.get_or_create(coefficient=seq.count('C')))
            reaction.participants.append(metabolites['gmp'][
                rna_compartment.id].species_coefficients.get_or_create(coefficient=seq.count('G')))
            reaction.participants.append(metabolites['ump'][
                rna_compartment.id].species_coefficients.get_or_create(coefficient=seq.count('U')))
            reaction.participants.append(metabolites['h'][
                rna_compartment.id].species_coefficients.get_or_create(coefficient=len(seq)-1))
                             
            # Assign modifier
            exo_ids = rna_exo_pair[rna_kb.id]
            modifier_obs = []
            for exo_id in exo_ids:
                complexes = model.species_types.get(name=exo_id)

                if not complexes:
                    raise ValueError('{} that catalyzes the degradation of {} cannot be found'.format(exo_id, rna_kb.id))

                else:                
                    observable = model.observables.get_one(
                        name='{} observable in {}'.format(exo_id, rna_compartment.name))

                    if not observable:
                        
                        all_species = {}                             
                        
                        for compl_variant in complexes:
                            exo_species = compl_variant.species.get_one(compartment=rna_compartment)
                            if not exo_species:
                                raise ValueError('{} cannot be found in the {}'.format(exo_species, rna_compartment.name))
                            all_species[exo_species.gen_id()] = exo_species
                            
                        observable_expression, error = wc_lang.ObservableExpression.deserialize(
                            ' + '.join(list(all_species.keys())), {
                            wc_lang.Species: all_species,
                            })
                        assert error is None, str(error)

                        observable = model.observables.create(
                                name='{} observable in {}'.format(exo_id, rna_compartment.name),
                                expression=observable_expression)
                        observable.id = 'obs_{}'.format(len(model.observables))

                if observable not in modifier_obs:
                    modifier_obs.append(observable)    

            if len(modifier_obs) == 1:
                self._degradation_modifier[reaction.name] = modifier_obs[0]
            else:
                all_obs = {obs.id: obs for obs in modifier_obs} 
                observable_expression, error = wc_lang.ObservableExpression.deserialize(
                    ' + '.join(list(all_obs.keys())), {
                    wc_lang.Observable: all_obs,
                    })
                assert error is None, str(error)

                observable = model.observables.create(
                            name='Combined exosome observable in {}'.format(rna_compartment.name),
                            expression=observable_expression)
                observable.id = 'obs_{}'.format(len(model.observables))

                self._degradation_modifier[reaction.name] = observable        
            
    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """
                
        for reaction in self.submodel.reactions:

            modifier = self._degradation_modifier[reaction.name]

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

        rnas_kb = cell.species_types.get(__type=wc_kb.eukaryote_schema.TranscriptSpeciesType)
        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):

            init_species_counts = {}
        
            modifier = self._degradation_modifier[reaction.name]      
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

            rna_reactant = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=rna_compartment)

            half_life = rna_kb.properties.get_one(property='half_life').get_value()
            mean_concentration = rna_reactant.distribution_init_concentration.mean

            average_rate = utils.calc_avg_deg_rate(mean_concentration, half_life)

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
