""" Generator for RNA degradation submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen
import wc_lang
import wc_kb
import numpy
import math

class RnaDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for RNA degradation submodel """

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='rna_degradation')
        cytosol = model.compartments.get_one(id='c')

        amp = model.species_types.get_one(id='amp').species.get_one(compartment=cytosol)
        cmp = model.species_types.get_one(id='cmp').species.get_one(compartment=cytosol)
        gmp = model.species_types.get_one(id='gmp').species.get_one(compartment=cytosol)
        ump = model.species_types.get_one(id='ump').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rna_kbs:

            rna_model = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            seq = rna_kb.get_seq()
            reaction = submodel.reactions.get_or_create(id=rna_kb.id.replace('rna_tu_', 'degrad_tu_'))
            reaction.name = rna_kb.id.replace('rna_', 'degrad_rna_')
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
            for degradosome_kb in cell.observables.get_one(id='degrade_rnase_obs').species:
                degradosome_species_type_model = model.species_types.get_one(id=degradosome_kb.species.species_type.id)
                degradosome_species_model = degradosome_species_type_model.species.get_one(compartment=cytosol)

                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(coefficient=(-1)*degradosome_kb.coefficient))
                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(coefficient=degradosome_kb.coefficient))

    def gen_phenom_rates(self):
        """ Generate rate laws with exponential dynamics """
        submodel = self.model.submodels.get_one(id='rna_degradation')
        cytosol = self.model.compartments.get_one(id='c')
        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        avg_rna_half_life = self.calc_mean_half_life(species_types_kb=rnas_kb)
        cell_cycle_length = self.knowledge_base.cell.properties.get_one(id='cell_cycle_length').value

        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):
            if (math.isnan(rna_kb.half_life) or rna_kb.half_life==0):
                half_life = avg_rna_half_life
            else:
                half_life = rna_kb.half_life

            specie_type_model = self.model.species_types.get_one(id=rna_kb.id)
            specie_model = specie_type_model.species.get_one(compartment=cytosol)

            rate_law = reaction.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = '({} / {}) * {}'.format(numpy.log(2), half_life, specie_model.id())

            rate_law.equation = wc_lang.RateLawEquation(expression = expression)
            rate_law.equation.modifiers.append(specie_model)

    def gen_mechanistic_rates(self):
        """ Generate rate laws with calibrated dynamics """
        submodel = self.model.submodels.get_one(id='rna_degradation')
        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        avg_rna_half_life = self.calc_mean_half_life(species_types_kb=rnas_kb)
        cell_cycle_length = self.knowledge_base.cell.properties.get_one(id='cell_cycle_length').value

        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):
            if (math.isnan(rna_kb.half_life) or rna_kb.half_life==0):
                half_life = avg_rna_half_life
            else:
                half_life = rna_kb.half_life

            self.gen_mechanistic_rate_law_eq(specie_type_kb=rna_kb,
                                             submodel=submodel,
                                             reaction=reaction,
                                             beta = 1,
                                             half_life=half_life,
                                             cell_cycle_length=cell_cycle_length)
