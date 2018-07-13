""" Generator for protein  degradation submodels based on KBs for random in silico organisms

:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
         Jonathan Karr <karr@mssm.edu>
:Date: 2018-07-05
:Copyright: 2018, Karr Lab
:License: MIT
"""
import numpy
import scipy
import wc_kb
import wc_lang
import wc_model_gen


class ProteinDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Gnerator for Protein degradation model"""

    def gen_compartments(self):
        self.cell = cell = self.knowledge_base.cell
        model = self.model
        self.cytosol = cytosol = model.compartments.get_or_create(id='c')[0]
        cytosol.name = 'cytosol'
        cytosol.initial_volume = cell.properties.get_one(
            id='mean_volume').value

    def gen_species(self):
        "Generate the protein species for the model"

        cell = self.cell
        model = self.model
        cytosol = self.cytosol

        proteins = cell.species_types.get(__type=wc_kb.ProteinSpeciesType)
        for protein in proteins:
            species_type = model.species_types.get_or_create(id=protein.id)
            if not species_type.name:
                species_type.name = protein.name
                species_type.type = wc_lang.SpeciesTypeType.protein
                species_type.structure = protein.get_seq()
                species_type.empirical_formula = protein.get_empirical_formula()
                species_type.molecular_weight = protein.get_mol_wt()
                species_type.charge = protein.get_charge()
                species = species_type.species.get_or_create(
                    compartment=cytosol)
                species.concentration = wc_lang.Concentration(
                    value=protein.concentration, units=wc_lang.ConcentrationUnit.M)
