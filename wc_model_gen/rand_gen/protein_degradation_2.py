""" Generator for protein  degradation submodels based on KBs for random in silico organisms

:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
         Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
         Balazs Szigeti <balazs.szigeti@mssm.edu>
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
    """ Generator for Protein degradation model"""

    def gen_compartments(self):
        self.cell = cell = self.knowledge_base.cell
        model = self.model
        self.cytosol = cytosol = model.compartments.get_or_create(id='c')
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

    def gen_reactions(self):
        submodel = self.submodel

        kb = self.knowledge_base

        cytosol = self.cytosol

        h2o = self.model.species_types.get_one(
            id='h2o').species.get_one(compartment=cytosol)
        h = self.model.species_types.get_one(
            id='h').species.get_one(compartment=cytosol)

        for specie_type in self.model.species_types.get(type=wc_lang.SpeciesTypeType.protein):

            if kb.cell.species_types.get_one(id=specie_type.id):
                # NoneType is species is not within KB
                # i.e. inactive species that represent some intermediary, e.g. _att (attached to ribosome)

                seq = kb.cell.species_types.get_one(id=specie_type.id).get_seq()

                reaction = wc_lang.core.Reaction(id='prot_degradation_' + specie_type.id, submodel=submodel)

                # Adding reaction participants LHS

                aaList = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", 
"M", "F", "P", "S", "T", "W", "Y", "V"]
                prot = specie_type.species.get_one(compartment=cytosol)
                reaction.participants = []
                reaction.participants.add(prot.species_coefficients.get_or_create(
                    coefficient=-1))
            
                for aa in aaList:
                    if aa in str(seq):
                        aa_species_type = self.model.species_types.get_or_create(id=aa)
                        aa_species = aa_species_type.species.get_or_create(compartment=cytosol)
                        reaction.participants.add(aa_species.species_coefficients.get_or_create(
                            coefficient=str(seq).count(aa)))
                        

                reaction.participants.add(h2o.species_coefficients.get_or_create(
                    coefficient=-(kb.cell.species_types.get_one(id=specie_type.id).get_len() - 1)))
                reaction.participants.add(h.species_coefficients.get_or_create(
                    coefficient=kb.cell.species_types.get_one(id=specie_type.id).get_len() - 1))

    def gen_rate_laws(self):
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = self.cytosol

        prots = cell.species_types.get(__type=wc_kb.ProteinSpeciesType)
        for prot, rxn in zip(prots, self.submodel.reactions):
            rl = rxn.rate_laws.create()
            rl.direction = wc_lang.RateLawDirection.forward
            rl.equation = wc_lang.RateLawEquation(
                expression='k_cat * {0}[c] / (k_m + {0}[c])'.format(prot.id))
            rl.k_cat = 1 #2 * numpy.log(2) / prot.half_life
            rl.k_m = 1 #prot.concentration
            rl.equation.modifiers.append(rxn.participants[0].species)

