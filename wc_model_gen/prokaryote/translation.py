""" Generating wc_lang formatted models from knowledge base.
:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen
import wc_lang
import wc_kb
import numpy
import math


class TranslationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generate translation submodel. """

    def gen_reactions(self):
        """ Generate a lumped reaction that cvers initiation, elongation and termination for each protein translated """
        model = self.model
        submodel = self.submodel
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        # Get species involved in reaction - tRna handeled on a per codon bases below
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        gdp = model.species_types.get_one(id='gdp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)
        initiation_factors = model.observables.get_one(id='translation_init_factors_obs').expression.species[0]
        elongation_factors = model.observables.get_one(id='translation_elongation_factors_obs').expression.species[0]
        release_factors = model.observables.get_one(id='translation_release_factors_obs').expression.species[0]

        bases = "TCAG"
        codons = [a + b + c for a in bases for b in bases for c in bases]

        proteins_kb = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb in proteins_kb:

            protein_model = model.species_types.get_one(id=protein_kb.id).species.get_one(compartment=cytosol)
            n_steps = protein_kb.get_len()
            reaction = model.reactions.get_or_create(submodel=submodel, id=protein_kb.id.replace('prot_', 'translation_'))
            reaction.name = protein_kb.id.replace('prot_', 'translation_')
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
                        trna = model.observables.get_one(id='tRNA_'+codon+'_obs').expression.species[0]
                        reaction.participants.add(trna.species_coefficients.get_or_create(coefficient=-n))

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
            for ribosome_kb in cell.observables.get_one(id='ribosome_obs').species:
                ribosome_species_type_model = model.species_types.get_one(id=ribosome_kb.species.species_type.id)
                ribosome_model = ribosome_species_type_model.species.get_one(compartment=cytosol)

                reaction.participants.add(ribosome_model.species_coefficients.get_or_create(coefficient=(-1)*ribosome_kb.coefficient))
                reaction.participants.add(ribosome_model.species_coefficients.get_or_create(coefficient=ribosome_kb.coefficient))

    def gen_phenom_rates(self):
        """ Generate rate laws with exponential dynamics """
        submodel = self.model.submodels.get_one(id='translation')
        proteins_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        cell_cycle_len = self.knowledge_base.cell.properties.get_one(id='cell_cycle_len').value

        for protein_kb, reaction in zip(proteins_kb, self.submodel.reactions):
            self.gen_phenom_rate_law_eq(specie_type_kb=protein_kb,
                                        reaction=reaction,
                                        half_life=protein_kb.half_life,
                                        cell_cycle_len=cell_cycle_len)

    def gen_mechanistic_rates(self):
        """ Generate rate laws associated with submodel """
        submodel = self.model.submodels.get_one(id='translation')
        proteins_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        cell_cycle_len = self.knowledge_base.cell.properties.get_one(id='cell_cycle_len').value

        for protein_kb, reaction in zip(proteins_kb, self.submodel.reactions):
            self.gen_mechanistic_rate_law_eq(specie_type_kb=protein_kb,
                                             submodel=submodel,
                                             reaction=reaction,
                                             beta=1,
                                             half_life=protein_kb.half_life,
                                             cell_cycle_len=cell_cycle_len)
