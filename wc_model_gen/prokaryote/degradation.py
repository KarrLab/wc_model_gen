""" Generating wc_lang formatted models from knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_lang
import wc_model_gen


class DegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):

    def gen_reactions(self):
        submodel = self.submodel
        species_type_types = [
            wc_lang.SpeciesTypeType.protein,
            wc_lang.SpeciesTypeType.rna,
        ]
        kb = self.knowledge_base

        for species_type_type in species_type_types:
            for specie_type in self.model.species_types.get(type=species_type_type):

                if kb.cell.species_types.get_one(id=specie_type.id):
                    # NoneType is species is not within KB
                    # i.e. inactive species that represent some intermediary, e.g. _att (attached to ribosome)

                    formula = kb.cell.species_types.get_one(id=specie_type.id).get_empirical_formula()

                    for specie in specie_type.species:
                        compartment = specie.compartment
                        reaction = wc_lang.core.Reaction(id='degradation_' + specie.species_type.id, submodel=submodel)

                        # Adding reaction participants LHS
                        reaction.participants.create(species=specie, coefficient=-1)

                        # Adding reaction participants RHS
                        for element in formula:
                            degrade_specie_type = self.model.species_types.get_or_create(id=element)
                            degrade_specie = degrade_specie_type.species.get_or_create(compartment=compartment)
                            reaction.participants.create(species=degrade_specie, coefficient=formula[element])

    def gen_rate_laws(self):
        submodel = self.submodel

        for reaction in submodel.reactions:
            exp = 'k_cat'
            mod = []

            for participant in reaction.participants:
                if participant.coefficient > 0:
                    continue

                if participant.coefficient < 0:
                    compartment = participant.species.compartment
                    exp = exp + ' * (' + participant.species.id() + '/ (k_m + ' + participant.species.id() + '))'
                    mod.append(self.model.species_types.get_one(
                        id=participant.species.species_type.id).species.get_one(compartment=compartment))

            rate_eq = wc_lang.core.RateLawEquation(expression=exp, modifiers=mod)
            rate_law = wc_lang.core.RateLaw(reaction=reaction,
                                            direction=wc_lang.core.RateLawDirection.forward,
                                            equation=rate_eq,
                                            k_cat=1,
                                            k_m=1)
