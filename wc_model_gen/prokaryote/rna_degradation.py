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


class RnaDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for RNA degradation submodel """

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

        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rna_kbs:

            rna_model = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            seq = rna_kb.get_seq()
            reaction = model.reactions.get_or_create(submodel=submodel, id=rna_kb.id.replace('rna_tu_', 'degrad_tu_'))
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

                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(
                    coefficient=(-1)*degradosome_kb.coefficient))
                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(
                    coefficient=degradosome_kb.coefficient))

    def gen_phenom_rates(self):
        """ Generate rate laws with exponential dynamics """
        model = self.model
        kb = self.knowledge_base
        submodel = model.submodels.get_one(id='rna_degradation')
        cytosol = model.compartments.get_one(id='c')
        rnas_kb = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        cell_cycle_len = kb.cell.properties.get_one(id='cell_cycle_len').value

        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):
            objects = {
                wc_lang.Species: {},
                wc_lang.Parameter: {},
            }

            species_type_model = model.species_types.get_one(id=rna_kb.id)
            species_model = species_type_model.species.get_one(compartment=cytosol)
            objects[wc_lang.Species][species_model.id] = species_model

            rate_law = model.rate_laws.create(
                id=wc_lang.RateLaw.gen_id(reaction.id, wc_lang.RateLawDirection.forward.name),
                reaction=reaction,
                direction=wc_lang.RateLawDirection.forward)

            half_life_model = model.parameters.get_or_create(id='half_life_{}'.format(species_type_model.id),
                                                             type=wc_lang.ParameterType.other,
                                                             value=rna_kb.half_life,
                                                             units='s')
            objects[wc_lang.Parameter][half_life_model.id] = half_life_model

            molecule_units = model.parameters.get_or_create(id='molecule_units',
                                                            type=wc_lang.ParameterType.other,
                                                            value=1.,
                                                            units='molecule')
            objects[wc_lang.Parameter][molecule_units.id] = molecule_units

            expression = '(log(2) / {}) / {} * {}'.format(half_life_model.id, molecule_units.id, species_model.id)
            rate_law.expression, error = wc_lang.RateLawExpression.deserialize(expression, objects)
            assert error is None, str(error)

    def gen_mechanistic_rates(self):
        """ Generate rate laws with calibrated dynamics """
        model = self.model
        kb = self.knowledge_base

        submodel = model.submodels.get_one(id='rna_degradation')
        rnas_kb = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        cell_cycle_len = kb.cell.properties.get_one(id='cell_cycle_len').value

        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):
            self.gen_mechanistic_rate_law_eq(specie_type_kb=rna_kb,
                                             submodel=submodel,
                                             reaction=reaction,
                                             beta=1.,
                                             half_life=rna_kb.half_life,
                                             cell_cycle_len=cell_cycle_len)
