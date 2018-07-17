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

        cell = self.knowledge_base.cell
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

    def gen_reactions(self):

        model = self.model
        submodel = self.submodel
        cytosol = model.compartments.get_one(id='c')

        h20 = model.species_types.get_one(
            id='h20').species_types.get_one(compartment=cytosol)
        atp = model.species_types.get_one(
            id='atp').species_types.get_one(compartment=cytosol)
        adp = model.species_types.get_one(
            id='adp').species_types.get_one(compartment=cytosol)
        pi = model.species_types.get_one(
            id='pi').species_types.get_one(compartment=cytosol)

        ala = model.species_types.get_one(
            id='ala').species_types.get_one(compartment=cytosol)
        arg = model.species_types.get_one(
            id='arg').species_types.get_one(compartment=cytosol)
        asp = model.species_types.get_one(
            id='asp').species_types.get_one(compartment=cytosol)
        asn = model.species_types.get_one(
            id='asn').species_types.get_one(compartment=cytosol)
        cys = model.species_types.get_one(
            id='cys').species_types.get_one(compartment=cytosol)
        gln = model.species_types.get_one(
            id='gln').species_types.get_one(compartment=cytosol)
        glu = model.species_types.get_one(
            id='glu').species_types.get_one(compartment=cytosol)
        gly = model.species_types.get_one(
            id='gly').species_types.get_one(compartment=cytosol)
        his = model.species_types.get_one(
            id='his').species_types.get_one(compartment=cytosol)
        ile = model.species_types.get_one(
            id='ile').species_types.get_one(compartment=cytosol)
        leu = model.species_types.get_one(
            id='leu').species_types.get_one(compartment=cytosol)
        lys = model.species_types.get_one(
            id='lys').species_types.get_one(compartment=cytosol)
        met = model.species_types.get_one(
            id='met').species_types.get_one(compartment=cytosol)
        phe = model.species_types.get_one(
            id='phe').species_types.get_one(compartment=cytosol)
        pro = model.species_types.get_one(
            id='pro').species_types.get_one(compartment=cytosol)
        ser = model.species_types.get_one(
            id='ser').species_types.get_one(compartment=cytosol)
        thr = model.species_types.get_one(
            id='thr').species_types.get_one(compartment=cytosol)
        trp = model.species_types.get_one(
            id='trp').species_types.get_one(compartment=cytosol)
        tyr = model.species_types.get_one(
            id='tyr').species_types.get_one(compartment=cytosol)
        val = model.species_types.get_one(
            id='val').species_types.get_one(compartment=cytosol)

        proteins = self.model.species_types.get(
            __type=wc_kb.core.ProteinSpeciesType)

        degradation_atpase = numpy.random.choice(proteins)
        while degradation_atpase.name:
            degradation_atpase = numpy.random.choice(proteins)

        degradation_protease = numpy.random.choice(proteins)
        while degradation_protease.name:
            degradation_protease = numpy.random.choice(proteins)

        for protein in proteins:
            rxn = submodel.reactions.get_or_create(
                id=protein.id.replace('protein_', 'protein_degradation_'))
            rxn.name = protein.id.replace('protein', 'protein degradation')
