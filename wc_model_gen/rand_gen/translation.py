""" Generating wc_lang formatted models from knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import wc_lang
import wc_model_gen
import random

class TranslationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generate translation submodel. """

    def gen_species(self):
        """ Generate a set of 2 protein species (_att: attached to ribosome and acitve form) for each protein """
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        # Initiate ribosome species types (complexes)
        species_type = self.model.species_types.create(id='complex_70S_IA', name='complex_70S_IA', type=5)
        species = species_type.species.create(compartment=compartment)
        species.concentration = wc_lang.core.Concentration(value=1e-2, units='M')

        species_type = self.model.species_types.create(id='complex_70S_A', name='complex_70S_A', type=5)
        species = species_type.species.create(compartment=compartment)
        species.concentration = wc_lang.core.Concentration(value=1e-2, units='M')

        # Create both functional and afunctional form (_att: attached to RNA) of every protein in KB
        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):

            # Add functional form of protein
            species_type = self.model.species_types.create(
                id=protein.id,
                type=2,
                name=protein.name,
                structure=protein.get_seq(),
                empirical_formula=str(protein.get_empirical_formula()),
                molecular_weight=str(protein.get_mol_wt()),
                charge=str(protein.get_charge()))

            species = species_type.species.create(compartment=compartment)
            species.concentration = wc_lang.core.Concentration(value=0, units='M')

            # Add inactive form of protein, attached to ribosome
            species_type = self.model.species_types.create(
                id=protein.id+'_att',
                type=2,
                name=protein.name+'_att',
                structure=protein.get_seq(),
                empirical_formula=str(protein.get_empirical_formula()),
                molecular_weight=str(protein.get_mol_wt()),
                charge=str(protein.get_charge()))

            species = species_type.species.create(compartment=compartment)
            species.concentration = wc_lang.core.Concentration(value=0, units='M')

    def gen_reactions(self):
        """ Generate a set of 3 reqactions (initation, elongation, termination) for each protein """

        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')
        prots = self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType)
        rnas = self.knowledge_base.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType)
        trnas = []
        for rna in rnas:
            if rna.type == wc_kb.RnaType.tRna:
                trnas.append(rna)
        IF1 = random.choice(prots).id
        IF2 = random.choice(prots).id
        IF3 = random.choice(prots).id
        EFtu = random.choice(prots).id
        EFts = random.choice(prots).id
        EFg = random.choice(prots).id
        RF = random.choice(prots).id
        amino_acids = {'S': ['tRNA-Ser', random.choice(trnas).id],
                       'L': ['tRNA-Leu', random.choice(trnas).id],
                       'R': ['tRNA-Arg', random.choice(trnas).id],
                       'T': ['tRNA-Thr', random.choice(trnas).id],
                       'G': ['tRNA-Gly', random.choice(trnas).id],
                       'F': ['tRNA-Phe', random.choice(trnas).id],
                       'W': ['tRNA-Trp', random.choice(trnas).id],
                       'K': ['tRNA-Lys', random.choice(trnas).id],
                       'I': ['tRNA-Ile', random.choice(trnas).id],
                       'A': ['tRNA-Ala', random.choice(trnas).id],
                       'M': ['tRNA-Met', random.choice(trnas).id],
                       'Q': ['tRNA-Gln', random.choice(trnas).id],
                       'E': ['tRNA-Glu', random.choice(trnas).id],
                       'P': ['tRNA-Pro', random.choice(trnas).id],
                       'V': ['tRNA-Val', random.choice(trnas).id],
                       'C': ['tRNA-Cys', random.choice(trnas).id],
                       'Y': ['tRNA-Tyr', random.choice(trnas).id],
                       'H': ['tRNA-His', random.choice(trnas).id],
                       'N': ['tRNA-Asn', random.choice(trnas).id],
                       'D': ['tRNA-Asp', random.choice(trnas).id]}
        # Add translation initating reactions
        for protein in prots:
            reaction = wc_lang.core.Reaction(id='translation_init_' + protein.id, submodel=submodel)

            # Adding reaction participants LHS
            specie = self.model.species_types.get_one(id='complex_70S_IA').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)
            specie = self.model.species_types.get_one(id='gtp').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)
            specie = self.model.species_types.get_one(id=IF1).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)
            specie = self.model.species_types.get_one(id=IF2).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)
            specie = self.model.species_types.get_one(id=IF3).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)

            # Adding reaction participants RHS
            specie = self.model.species_types.get_one(id='complex_70S_A').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)
            specie = self.model.species_types.get_one(id='gdp').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)
            specie = self.model.species_types.get_one(id='pi').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)
            specie = self.model.species_types.get_one(id=IF1).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)
            specie = self.model.species_types.get_one(id=IF2).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)
            specie = self.model.species_types.get_one(id=IF3).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)

        # Add translation elongation reactions
        for protein in prots:
            reaction = wc_lang.core.Reaction(id='translation_elon_' + protein.id, submodel=submodel)
            n_steps = len(protein.get_seq())

            # Adding reaction participants LHS
            specie = self.model.species_types.get_one(id='complex_70S_A').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)
            specie = self.model.species_types.get_one(id=EFtu).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-n_steps)
            specie = self.model.species_types.get_one(id=EFts).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-n_steps)
            specie = self.model.species_types.get_one(id=EFg).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-n_steps)
            specie = self.model.species_types.get_one(id='gtp').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-2*n_steps)

            # tRNAs - add more testing to this bit
                    
            for letter in list(set(str(protein.get_seq()))):
                if letter == '*':
                    continue

                n = protein.get_seq().tostring().count(letter)
                specie = self.model.species_types.get_or_create(
                    id=amino_acids[letter][1], type=4).species.get_or_create(compartment=compartment)
                reaction.participants.create(species=specie, coefficient=-n)

            # Adding reaction participants RHS
            specie = self.model.species_types.get_one(id=protein.id+'_att').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)
            specie = self.model.species_types.get_one(id='gdp').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=2*n_steps)
            specie = self.model.species_types.get_one(id='pi').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=2*n_steps)

        # Add translation termination reactions
        for protein in prots:
            reaction = wc_lang.core.Reaction(id='translation_term_' + protein.id, submodel=submodel)
            n_steps = len(protein.get_seq())

            # Adding reaction participants LHS
            specie = self.model.species_types.get_one(id=RF).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)
            specie = self.model.species_types.get_one(id='gtp').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)
            specie = self.model.species_types.get_one(id=protein.id+'_att').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=-1)

            # Adding reaction participants RHS
            specie = self.model.species_types.get_one(id=protein.id).species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)
            specie = self.model.species_types.get_one(id='gdp').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)
            specie = self.model.species_types.get_one(id='pi').species.get_one(compartment=compartment)
            reaction.participants.create(species=specie, coefficient=1)

    def gen_rate_laws(self):
        """ Generate the rate laws associate dwith reactions """
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        for reaction in submodel.reactions:
            exp = 'k_cat'
            mod = []

            for participant in reaction.participants:
                if participant.coefficient > 0:
                    continue

                if participant.coefficient < 0:
                    exp = exp + ' * (' + participant.species.id() + '/ (k_m + ' + participant.species.id() + '))'
                    mod.append(self.model.species_types.get_one(
                        id=participant.species.species_type.id).species.get_one(compartment=compartment))

            rate_eq = wc_lang.core.RateLawEquation(expression=exp, modifiers=mod)
            rate_law = wc_lang.core.RateLaw(reaction=reaction,
                                            direction=wc_lang.core.RateLawDirection.forward,
                                            equation=rate_eq,
                                            k_cat=1,
                                            k_m=1)
