""" Generator for RNA degradation submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import wc_model_gen.utils as utils
import scipy.constants
import wc_model_gen
import wc_lang
import wc_kb
import numpy
import math


class RnaDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for RNA degradation submodel

        Options:
        * beta (:obj:`float`, optional): ratio of Michaelis-Menten constant to substrate
            concentration (Km/[S]) for use when estimating Km values, the default value is 1
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        beta = options.get('beta', 1.)
        options['beta'] = beta

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        kb = self.knowledge_base
        cell = kb.cell
        submodel = model.submodels.get_one(id='rna_degradation')
        cytosol = model.compartments.get_one(id='c')

        amp = model.species_types.get_one(id='amp').species.get_one(compartment=cytosol)
        cmp = model.species_types.get_one(id='cmp').species.get_one(compartment=cytosol)
        gmp = model.species_types.get_one(id='gmp').species.get_one(compartment=cytosol)
        ump = model.species_types.get_one(id='ump').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote.RnaSpeciesType)
        for rna_kb in rna_kbs:

            rna_model = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            seq = rna_kb.get_seq()
            reaction = model.reactions.get_or_create(submodel=submodel, id='degradation_' + rna_kb.id)
            reaction.name = 'degradation ' + rna_kb.name
            reaction.participants = []

            # Adding participants to LHS
            reaction.participants.add(rna_model.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.add(h2o.species_coefficients.get_or_create(coefficient=-(rna_kb.get_len() - 1)))

            # Adding participants to RHS
            reaction.participants.add(amp.species_coefficients.get_or_create(coefficient=seq.count('A')))
            reaction.participants.add(cmp.species_coefficients.get_or_create(coefficient=seq.count('C')))
            reaction.participants.add(gmp.species_coefficients.get_or_create(coefficient=seq.count('G')))
            reaction.participants.add(ump.species_coefficients.get_or_create(coefficient=seq.count('U')))
            reaction.participants.add(h.species_coefficients.get_or_create(coefficient=rna_kb.get_len() - 1))

            # Add members of the degradosome
            # Counterintuitively .specie is a KB species_coefficient object
            for degradosome_kb in cell.observables.get_one(id='degrade_rnase_obs').expression.species:
                degradosome_species_type_model = model.species_types.get_one(id=degradosome_kb.species_type.id)
                degradosome_species_model = degradosome_species_type_model.species.get_one(compartment=cytosol)

                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(
                    coefficient=-1))
                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(
                    coefficient=1))

    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """

        model = self.model

        modifier = model.observables.get_one(id='degrade_rnase_obs')

        for reaction in self.submodel.reactions:

            rate_law_exp, parameters = utils.gen_michaelis_menten_like_rate_law(
                model, reaction, modifiers=[modifier])
            
            rate_law = model.rate_laws.create(direction=wc_lang.RateLawDirection.forward,
                                              type=None,
                                              expression=rate_law_exp,
                                              reaction=reaction,
                                              )
            rate_law.id = rate_law.gen_id()

    def calibrate_submodel(self):
        """ Calibrate the submodel using data in the KB """
        model = self.model
        beta = self.options.get('beta')

        Avogadro = model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))

        cytosol = model.compartments.get_one(id='c')

        init_species_counts = {}

        modifier = model.observables.get_one(id='degrade_rnase_obs')
        for species in modifier.expression.species:
            init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean

        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote.RnaSpeciesType)
        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):

            rna_reactant = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            half_life = rna_kb.properties.get_one(property='half_life').get_value()
            mean_concentration = rna_reactant.distribution_init_concentration.mean

            average_rate = average_rate = utils.calc_avg_deg_rate(mean_concentration, half_life)

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
                wc_lang.Compartment: {cytosol.id: cytosol.init_volume.mean * cytosol.init_density.value},
            })
