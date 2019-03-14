""" Generator for protein  degradation submodels based on KBs for random in silico organisms

:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-07-05
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_utils.util.ontology import wcm_ontology
from wc_utils.util.units import unit_registry
import wc_model_gen.utils as utils
import scipy.constants
import wc_model_gen
import wc_lang
import wc_kb
import numpy
import math


class ProteinDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for protein degradation model 

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
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='protein_degradation')
        cytosol = model.compartments.get_one(id='c')

        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        adp = model.species_types.get_one(id='adp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)

        amino_acids = ['ala', 'arg', 'asp', 'asn', 'cys', 'gln', 'glu', 'gly', 'his',
                       'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val']

        aas = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
               "S", "T", "W", "Y", "V"]

        proteins_kb = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)

        for protein_kb in proteins_kb:

            protein_model = model.species_types.get_one(id=protein_kb.id).species.get_one(compartment=cytosol)
            seq = protein_kb.get_seq()
            reaction = model.reactions.get_or_create(submodel=submodel, id='degradation_' + protein_kb.id)
            reaction.name = 'degradation ' + protein_kb.name
            reaction.participants = []

            # Adding participants to LHS
            reaction.participants.add(protein_model.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.add(atp.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.add(h2o.species_coefficients.get_or_create(coefficient=-(len(seq)-1)))

            # Adding participants to RHS
            reaction.participants.add(adp.species_coefficients.get_or_create(coefficient=1))
            reaction.participants.add(pi.species_coefficients.get_or_create(coefficient=1))

            # The code below should be used as currently tRNAs and AAs are always associated
            codons = []
            for start_position in range(0, len(protein_kb.gene.get_seq())-3, 3):
                codons.append(str(protein_kb.gene.get_seq()[start_position:start_position+3]))

            for codon in set(codons):
                obs_model_id = 'tRNA_' + codon + '_obs'
                obs_model = model.observables.get_one(id=obs_model_id)
                for specie in obs_model.expression.species:
                    reaction.participants.add(
                        specie.species_coefficients.get_or_create(coefficient=codons.count(codon)))

            # for amino_acid, aa in zip(amino_acids, aas):
            #    species = model.species_types.get_one(id=amino_acid).species.get_one(compartment=cytosol)
            #    reaction.participants.add(species.species_coefficients.get_or_create(coefficient=seq.count(aa)))

            # Add members of the degradosome
            # Counterintuitively .specie is a KB species_coefficient object
            for degradosome_kb in cell.observables.get_one(id='degrade_protease_obs').expression.species:
                degradosome_species_type_model = model.species_types.get_one(id=degradosome_kb.species_type.id)
                degradosome_species_model = degradosome_species_type_model.species.get_one(compartment=cytosol)

                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(
                    coefficient=-1))
                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(
                    coefficient=1))

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

        modifier = model.observables.get_one(id='degrade_protease_obs')

        for reaction in self.submodel.reactions:
            modifier_reactant = [i for i in modifier.expression.species if i.species_type.id in reaction.id]
            if modifier_reactant:
                rate_law_exp, parameters = utils.gen_michaelis_menten_like_rate_law(
                    Avogadro, molecule_units, reaction, modifiers=[modifier],
                    modifier_reactants=modifier_reactant)
            else:
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

        init_species_counts = {}

        modifier = model.observables.get_one(id='degrade_protease_obs')
        for species in modifier.expression.species:
            init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean

        proteins_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb, reaction in zip(proteins_kb, self.submodel.reactions):

            protein_reactant = model.species_types.get_one(id=protein_kb.id).species.get_one(compartment=cytosol)
            half_life = protein_kb.half_life
            mean_concentration = protein_reactant.distribution_init_concentration.mean

            average_rate = utils.calc_avg_deg_rate(mean_concentration, half_life)

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
