""" Generating wc_lang formatted models from knowledge base.
:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-01-21
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


class TranslationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generate translation submodel 

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
        """ Generate a lumped reaction that covers initiation, elongation and termination for each protein translated """
        model = self.model
        submodel = self.submodel
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        # Get species involved in reaction - tRna handeled on a per codon bases below
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        gdp = model.species_types.get_one(id='gdp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)

        # Get initiation factors, elongation factors, and release factors (modifiers)
        self._modifiers = [model.observables.get_one(id='translation_init_factors_obs'),
                           model.observables.get_one(id='translation_elongation_factors_obs'),
                           model.observables.get_one(id='translation_release_factors_obs')]
        initiation_factors = self._modifiers[0].expression.species[0]
        elongation_factors = self._modifiers[1].expression.species[0]
        release_factors = self._modifiers[2].expression.species[0]

        # Check if initiation, elongation & termination factors consist of one specie
        for modifier in self._modifiers:
            assert(len(modifier.expression.species) == 1)

        bases = "TCAG"
        codons = [a + b + c for a in bases for b in bases for c in bases]

        proteins_kb = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb in proteins_kb:

            protein_model = model.species_types.get_one(id=protein_kb.id).species.get_one(compartment=cytosol)
            n_steps = protein_kb.get_len()
            reaction = model.reactions.get_or_create(submodel=submodel, id='translation_' + protein_kb.id)
            reaction.name = 'translation ' + protein_kb.name
            reaction.participants = []

            # Adding participants to LHS
            reaction.participants.add(gtp.species_coefficients.get_or_create(coefficient=-(n_steps+2)))
            reaction.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=-n_steps))
            reaction.participants.add(release_factors.species_coefficients.get_or_create(coefficient=-1))

            # Add tRNAs to LHS
            for codon in codons:
                if codon not in ['TAG', 'TAA', 'TGA']:
                    n = 0

                    for base in range(0, len(protein_kb.gene.get_seq()), 3):
                        n += str(protein_kb.gene.get_seq()[base:base+3]).count(codon)

                    if n > 0:
                        trna_obs = model.observables.get_one(id='tRNA_'+codon+'_obs')
                        if trna_obs not in self._modifiers:
                            self._modifiers.append(trna_obs)
                        trna = trna_obs.expression.species.get_one(compartment=cytosol)

                        # tRNAs are modifiers
                        reaction.participants.add(trna.species_coefficients.get_or_create(coefficient=-n))
                        reaction.participants.add(trna.species_coefficients.get_or_create(coefficient=n))

                        # Add appropiate amino acids
                        # TODO: add in all AAs
                        if codon == 'ATG':
                            aa = model.species_types.get_one(id='met').species.get_one(compartment=cytosol)
                        elif codon == 'ACT' or codon == 'ACC' or codon == 'ACA' or codon == 'ACG':
                            aa = model.species_types.get_one(id='thr').species.get_one(compartment=cytosol)
                        elif codon == 'ATT' or codon == 'ATC' or codon == 'ATA':
                            aa = model.species_types.get_one(id='ile').species.get_one(compartment=cytosol)
                        elif codon == 'TTA' or codon == 'TTG' or codon == 'CTT' or codon == 'CTC' or codon == 'CTA' or codon == 'CTG':
                            aa = model.species_types.get_one(id='leu').species.get_one(compartment=cytosol)
                        else:
                            raise ValueError('Unknown codon: {}'.format(codon))

                        reaction.participants.add(aa.species_coefficients.get_or_create(coefficient=-n))

            # Adding participants to RHS
            if protein_model == initiation_factors:
                reaction.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=2))
                reaction.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps))
                reaction.participants.add(release_factors.species_coefficients.get_or_create(coefficient=1))

            elif protein_model == elongation_factors:
                reaction.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps+1))
                reaction.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=1))
                reaction.participants.add(release_factors.species_coefficients.get_or_create(coefficient=1))

            elif protein_model == release_factors:
                reaction.participants.add(release_factors.species_coefficients.get_or_create(coefficient=2))
                reaction.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=1))
                reaction.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps))

            else:
                reaction.participants.add(protein_model.species_coefficients.get_or_create(coefficient=1))
                reaction.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=1))
                reaction.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps))
                reaction.participants.add(release_factors.species_coefficients.get_or_create(coefficient=1))

            reaction.participants.add(gdp.species_coefficients.get_or_create(coefficient=n_steps+2))
            reaction.participants.add(pi.species_coefficients.get_or_create(coefficient=2*n_steps))

            # Add ribosome
            if model.observables.get_one(id='ribosome_obs') not in self._modifiers:
                self._modifiers.append(model.observables.get_one(id='ribosome_obs'))

            for ribosome_kb in cell.observables.get_one(id='ribosome_obs').expression.species:
                ribosome_species_type_model = model.species_types.get_one(id=ribosome_kb.species_type.id)
                ribosome_model = ribosome_species_type_model.species.get_one(compartment=cytosol)

                reaction.participants.add(ribosome_model.species_coefficients.get_or_create(coefficient=-1))
                reaction.participants.add(ribosome_model.species_coefficients.get_or_create(coefficient=1))

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

        for reaction in self.submodel.reactions:
            rate_law_exp, parameters = utils.gen_michaelis_menten_like_rate_law(
                Avogadro, molecule_units, reaction, modifiers=self._modifiers)
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

        for modifier in self._modifiers:
            for species in modifier.expression.species:
                init_species_counts[species.gen_id()] = species.distribution_init_concentration.mean

        protein_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb, reaction in zip(protein_kb, self.submodel.reactions):

            protein_product = model.species_types.get_one(id=protein_kb.id).species.get_one(compartment=cytosol)
            half_life = protein_kb.half_life
            mean_concentration = protein_product.distribution_init_concentration.mean

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
