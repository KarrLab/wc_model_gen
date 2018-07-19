""" Generator for metabolism submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen
import wc_lang


class MetabolismSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for metabolism submodel """

    def gen_compartments(self):
        cell = self.knowledge_base.cell
        model = self.model

        # Generate the compartments that are defined in the knowledge base
        for kb_comp in cell.compartments:
            model_comp = model.compartments.get_or_create(kb_comp.id)
            model_comp.name = kb_comp.name

        # If the kb defines "c" and "e" compartments, their properties will be set. If not, they will be created
        cyt = model.compartments.get_or_create(id='c')
        if not cyt.name:
            cyt.name = 'cytosol'
        cyt.initial_volume = cell.properties.get_one(id='mean_volume').value

        ext = model.compartments.get_or_create(id='e')
        if not ext.name:
            ext.name = 'extracellular space'
        ext.initial_volume = 1. / \
            cell.properties.get_one(id='mean_cell_density').value

    def gen_parameters(self):
        cell = self.knowledge_base.cell
        model = self.model
        param = model.parameters.get_or_create(id='fractionDryWeight')
        param.submodels.append(self.submodel)
        param.value = cell.properties.get_one(
            id='mean_fraction_dry_weight').value
        param.units = 'dimensionless'

    def gen_species(self):
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

        # get or create metabolite species
        ids = ['atp', 'ctp', 'gtp', 'utp', 'ppi', 'h2o',
               'h', 'amp', 'adp', 'cmp', 'gmp', 'ump', 'pi', 'gdp', 'ala', 'arg', 'asp', 'asn', 'cys', 'gln', 'glu', 'gly', 'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val']

        for id in ids:
            kb_met = cell.species_types.get_one(id=id)
            species_type = model.species_types.get_or_create(id=id)
            if not species_type.name:
                species_type.name = kb_met.id
                species_type.type = wc_lang.SpeciesTypeType.metabolite
                species_type.structure = kb_met.structure
                species_type.empirical_formula = kb_met.get_empirical_formula()
                species_type.molecular_weight = kb_met.get_mol_wt()
                species_type.charge = kb_met.get_charge()
                species_type_c = species_type.species.get_or_create(
                    compartment=cytosol)
                species_type_c.concentration = wc_lang.Concentration(
                    value=kb_met.concentration, units=wc_lang.ConcentrationUnit.M)

    def gen_reactions(self):
        pass

    def gen_rate_laws(self):
        pass
