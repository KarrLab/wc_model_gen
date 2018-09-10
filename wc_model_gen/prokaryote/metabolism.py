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

    def gen_a_specie(self, kb_metabolite, lang_compartment):
        """ Generate a species in a particular compartment

        Args:
            kb_metabolite (:obj:`wc_kb.MetaboliteSpeciesType`): a knowledgebase metabolite
            lang_compartment (:obj:`wc_lang.Compartment`): the wc_lang compartment containing the species

        Returns:
            :obj:`wc_lang.Species`: the species that was found or created
        """
        species_type = self.model.species_types.get_or_create(id=kb_metabolite.id)
        species_type.name = kb_metabolite.name
        species_type.type = wc_lang.SpeciesTypeType.metabolite
        species_type.structure = kb_metabolite.structure
        species_type.empirical_formula = kb_metabolite.get_empirical_formula()
        species_type.molecular_weight = kb_metabolite.get_mol_wt()
        species_type.charge = kb_metabolite.get_charge()
        species_type.comments = kb_metabolite.comments
        species = species_type.species.get_or_create(
            compartment=lang_compartment)
        species.concentration = wc_lang.Concentration(
            value=kb_metabolite.concentration, units=wc_lang.ConcentrationUnit.M)
        return species

    def gen_species(self):
        """ Generate all metabolic species in the cytosol """
        cytosol = self.model.compartments.get_one(id='c')
        metabolites = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.core.MetaboliteSpeciesType)

        # get or create metabolite species
        for kb_met in metabolites:
            self.gen_a_specie(kb_met, cytosol)

    def get_species_type_types(self, kb_rxn):
        """ Obtain the species type types used by a kb reaction """
        species_type_types = set()
        for participant in kb_rxn.participants:
            kb_species_type = participant.species.species_type
            species_type_types.add(type(kb_species_type))
        return species_type_types

    def gen_reactions(self):
        """ Generate reactions encoded within KB """

        for kb_rxn in self.knowledge_base.cell.reactions:
            # if species are metabolites, create lang reaction in metabolism
            # todo: generalize to all submodels
            if self.get_species_type_types(kb_rxn) == set([wc_kb.core.MetaboliteSpeciesType]):
                lang_rxn = self.submodel.reactions.create(
                    id=kb_rxn.id,
                    name=kb_rxn.name,
                    reversible=kb_rxn.reversible,
                    comments=kb_rxn.comments)
                for participant in kb_rxn.participants:
                    kb_species = participant.species
                    lang_species_type = self.model.species_types.get_one(
                        id=kb_species.species_type.id)
                    lang_compartment = self.model.compartments.get_one(
                        id=kb_species.compartment.id)
                    lang_species = lang_species_type.species.get_one(
                        compartment=lang_compartment)
                    # ensure that species are present in extracellular space
                    if lang_species is None:
                        lang_species = self.gen_a_specie(kb_species.species_type, lang_compartment)
                    lang_rxn.participants.add(
                        lang_species.species_coefficients.get_or_create(
                            coefficient=participant.coefficient))

    def gen_rate_laws(self):

        """ Generate rate laws for reactions encoded in KB """
        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')

        for kb_rxn in self.knowledge_base.cell.reactions:
            if self.get_species_type_types(kb_rxn) == set([wc_kb.core.MetaboliteSpeciesType]):
                lang_rxn = self.submodel.reactions.get_one(id=kb_rxn.id)
                for kb_rate_law in kb_rxn.rate_laws:
                    lang_rate_law = wc_lang.RateLaw(k_cat=kb_rate_law.k_cat,
                        k_m=kb_rate_law.k_m,
                        comments=kb_rate_law.comments)
                    lang_rxn.rate_laws.add(lang_rate_law)
                    if kb_rate_law.direction == wc_kb.RateLawDirection.forward:
                        lang_rate_law.direction = wc_lang.RateLawDirection.forward
                    elif kb_rate_law.direction == wc_kb.RateLawDirection.backward:  # pragma branch not covered
                        lang_rate_law.direction = wc_lang.RateLawDirection.backward
                    lang_rate_law.equation = wc_lang.RateLawEquation(
                        expression=kb_rate_law.equation.expression
                    )

                    for kb_modifier in kb_rate_law.equation.modifiers:
                        lang_species_type = self.model.species_types.get_one(
                            id=kb_modifier.species_type.id)
                        lang_species = lang_species_type.species.get_one(
                            compartment=c)
                        lang_rate_law.equation.modifiers.add(lang_species)
