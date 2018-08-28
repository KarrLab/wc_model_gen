""" Generator for metabolism submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen
import wc_lang
import wc_kb


class MetabolismSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for metabolism submodel """

    def gen_compartments(self):
        cell = self.knowledge_base.cell
        model = self.model

        # Generate the compartments that are defined in the knowledge base
        for kb_comp in cell.compartments:
            model_comp = model.compartments.get_or_create(id=kb_comp.id)
            model_comp.name = kb_comp.name

        # If the kb defines "c" and "e" compartments, their properties will be set. If not, they will be created
        cyt = model.compartments.get_or_create(id='c')
        if not cyt.name:
            cyt.name = 'cytosol'
        cyt.initial_volume = cell.properties.get_one(id='initial_volume').value

        ext = model.compartments.get_or_create(id='e')
        if not ext.name:
            ext.name = 'extracellular space'
        ext.initial_volume = 1. / cell.properties.get_one(id='cell_density').value

    def gen_parameters(self):
        cell = self.knowledge_base.cell
        model = self.model
        param = model.parameters.get_or_create(id='fraction_dry_weight')
        param.submodels.append(self.submodel)
        param.value = cell.properties.get_one(
            id='fraction_dry_weight').value
        param.units = 'dimensionless'

    def gen_species(self):
        """ Generate species associated with submodel """
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]
        metabolites = cell.species_types.get(
            __type=wc_kb.core.MetaboliteSpeciesType)

        # get or create metabolite species
        for kb_met in metabolites:
            species_type = model.species_types.get_or_create(id=kb_met.id)
            if not species_type.name:
                species_type.name = kb_met.name
                species_type.type = wc_lang.SpeciesTypeType.metabolite
                species_type.structure = kb_met.structure
                species_type.empirical_formula = kb_met.get_empirical_formula()
                species_type.molecular_weight = kb_met.get_mol_wt()
                species_type.charge = kb_met.get_charge()
                species = species_type.species.get_or_create(
                    compartment=cytosol)
                species.concentration = wc_lang.Concentration(
                    value=kb_met.concentration, units=wc_lang.ConcentrationUnit.M)

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        cytosol = self.model.compartments.get_one(id='c')

        for kb_rxn in self.knowledge_base.cell.reactions:
            lang_rxn = self.submodel.reactions.create(
                id=kb_rxn.id,
                name=kb_rxn.name,
                reversible=kb_rxn.reversible,
                comments=kb_rxn.comments)
            for participant in kb_rxn.participants:
                lang_species_type = self.model.species_types.get_one(
                    id=participant.species.species_type.id)
                lang_species = lang_species_type.species.get_one(
                    compartment=cytosol)
                lang_rxn.participants.add(
                    lang_species.species_coefficients.get_or_create(
                        coefficient=participant.coefficient))

    def gen_rate_laws(self):
        """ Generate rate laws for reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        for kb_rxn in self.knowledge_base.cell.reactions:
            lang_rxn = self.submodel.reactions.get_one(id=kb_rxn.id)
            for kb_rate_law in kb_rxn.rate_laws:
                lang_rate_law = wc_lang.RateLaw(k_cat=kb_rate_law.k_cat,
                    k_m=kb_rate_law.k_m,
                    comments=kb_rate_law.comments)
                lang_rxn.rate_laws.add(lang_rate_law)
                if kb_rate_law.direction == wc_kb.RateLawDirection.forward:
                    lang_rate_law.direction = wc_lang.RateLawDirection.forward
                elif kb_rate_law.direction == wc_kb.RateLawDirection.backward:
                    lang_rate_law.direction = wc_lang.RateLawDirection.backward
                lang_rate_law.equation = wc_lang.RateLawEquation(
                    expression=kb_rate_law.equation.expression
                )
                for kb_modifier in kb_rate_law.equation.modifiers:
                    lang_species_type = self.model.species_types.get_one(
                        id=kb_modifier.species_type.id)
                    lang_species = lang_species_type.species.get_one(
                        compartment=cytosol)
                    lang_rate_law.equation.modifiers.add(lang_species)
