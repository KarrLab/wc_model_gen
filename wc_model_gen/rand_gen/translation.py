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
        species_type.molecular_weight = '1' #placeholder
        species = species_type.species.create(compartment=compartment)
        species.concentration = wc_lang.core.Concentration(value=1e-2, units='M')

        species_type = self.model.species_types.create(id='complex_70S_A', name='complex_70S_A', type=5)
        species_type.molecular_weight = '1' #placeholder
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
        protsList = []
        IF1 = random.choice(prots).id
        protsList.append(IF1)
        IF2 = random.choice(prots).id
        while IF2 in protsList:
            IF2 = random.choice(prots).id
        protsList.append(IF2)
        IF3 = random.choice(prots).id
        while IF3 in protsList:
            IF3 = random.choice(prots).id
        protsList.append(IF3)
        EFtu = random.choice(prots).id
        while EFtu in protsList:
            EFtu = random.choice(prots).id
        protsList.append(EFtu)
        EFts = random.choice(prots).id
        while EFts in protsList:
            EFts = random.choice(prots).id
        protsList.append(EFts)
        EFg = random.choice(prots).id
        while EFg in protsList:
            EFg = random.choice(prots).id
        protsList.append(EFg)
        RF = random.choice(prots).id
        while RF in protsList:
            IF2 = random.choice(prots).id
        protsList.append(IF2)
        amino_acids = {}
        trna = random.choice(trnas)
        amino_acids['S'] = ['tRNA-Ser', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['L'] = ['tRNA-Leu', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['R'] = ['tRNA-Arg', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['T'] = ['tRNA-Thr', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['G'] = ['tRNA-Gly', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['F'] = ['tRNA-Phe', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['W'] = ['tRNA-Trp', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['K'] = ['tRNA-Lys', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['I'] = ['tRNA-Ile', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['A'] = ['tRNA-Ala', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['M'] = ['tRNA-Met', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['Q'] = ['tRNA-Gln', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['E'] = ['tRNA-Glu', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['P'] = ['tRNA-Pro', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['V'] = ['tRNA-Val', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['C'] = ['tRNA-Cys', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['Y'] = ['tRNA-Tyr', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['H'] = ['tRNA-His', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['N'] = ['tRNA-Asn', trna.id]
        trnas.remove(trna)
        trna = random.choice(trnas)
        amino_acids['D'] = ['tRNA-Asp', trna.id]
        trnas.remove(trna)
        
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
