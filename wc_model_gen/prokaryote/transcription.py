""" Generator for transcription submodels based on KBs

:Author: Jonathan Karr <karr@mssm.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_utils.util.ontology import wcm_ontology
from wc_utils.util.units import unit_registry
import wc_model_gen.utils as utils
import math
import numpy
import scipy.constants
import wc_model_gen
import wc_lang
import wc_kb


class TranscriptionSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for transcription submodel 

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
        submodel = self.submodel
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        # Get species involved in reaction
        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        # Create reaction for each RNA
        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rna_kbs:

            rna_model = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            reaction = model.reactions.get_or_create(submodel=submodel, id='transcription_' + rna_kb.id)
            reaction.name = 'transcription ' + rna_kb.name
            reaction.participants = []
            seq = rna_kb.get_seq()

            # Adding participants to LHS
            reaction.participants.add(atp.species_coefficients.get_or_create(coefficient=-seq.count('A')))
            reaction.participants.add(ctp.species_coefficients.get_or_create(coefficient=-seq.count('C')))
            reaction.participants.add(gtp.species_coefficients.get_or_create(coefficient=-seq.count('G')))
            reaction.participants.add(utp.species_coefficients.get_or_create(coefficient=-seq.count('U')))
            # reaction.participants.add(h2o.species_coefficients.get_or_create(coefficient=-1))

            # Adding participants to RHS
            reaction.participants.add(rna_model.species_coefficients.get_or_create(coefficient=1))
            reaction.participants.add(ppi.species_coefficients.get_or_create(coefficient=rna_kb.get_len()-1))
            #reaction.participants.add(h.species_coefficients.get_or_create(coefficient=1 + rna_kb.get_len()))

            # Add RNA polymerease
            for rnap_kb in cell.observables.get_one(id='rna_polymerase_obs').expression.species:
                rnap_species_type_model = model.species_types.get_one(id=rnap_kb.species_type.id)
                rnap_model = rnap_species_type_model.species.get_one(compartment=cytosol)

                reaction.participants.add(rnap_model.species_coefficients.get_or_create(coefficient=-1))
                reaction.participants.add(rnap_model.species_coefficients.get_or_create(coefficient=1))

    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """
        model = self.model

        Avogadro = model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))
        molecule_units = model.parameters.get_or_create(
            id='molecule_units',
            type=None,
            value=1.,
            units=unit_registry.parse_units('molecule'))

        modifier = model.observables.get_one(id='rna_polymerase_obs')

        for reaction in self.submodel.reactions:
            rate_law_exp, parameters = utils.gen_michaelis_menten_like_rate_law(
                Avogadro, molecule_units, reaction, modifiers=[modifier])
            model.parameters += parameters

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

        mean_doubling_time = self.knowledge_base.cell.properties.get_one(id='mean_doubling_time').value

        init_species_counts = {}

        modifier = model.observables.get_one(id='rna_polymerase_obs')
        for species in modifier.expression.species:
            init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean

        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):

            rna_product = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            half_life = rna_kb.half_life
            mean_concentration = rna_product.distribution_init_concentration.mean

            average_rate = utils.calc_avg_syn_rate(
                mean_concentration, half_life, mean_doubling_time)

            for species in reaction.get_reactants():

                init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean

                if model.parameters.get(id='K_m_{}_{}'.format(reaction.id, species.species_type.id)):
                    model_Km = model.parameters.get_one(
                        id='K_m_{}_{}'.format(reaction.id, species.species_type.id))
                    model_Km.value = beta * species.distribution_init_concentration.mean \
                        / Avogadro.value / species.compartment.mean_init_volume

            model_kcat = model.parameters.get_one(id='k_cat_{}'.format(reaction.id))
            model_kcat.value = 1.
            model_kcat.value = average_rate / reaction.rate_laws[0].expression._parsed_expression.eval({
                wc_lang.Species: init_species_counts,
                wc_lang.Compartment: {cytosol.id: cytosol.mean_init_volume * cytosol.init_density.value},
            })
