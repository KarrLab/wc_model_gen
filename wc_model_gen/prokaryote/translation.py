""" Generating wc_lang formatted models from knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
- protein.get_seq() is throwing error, but protein.rna.get_seq().translate() works, look into

"""

import wc_kb
import wc_lang
import wc_model_gen


class TranslationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generate translation submodel. """

    def gen_species(self):
        """ Generate a set of 2 protein species (_att: attached to ribosome and acitve form) for each protein """
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        # Initiate ribosome species types (complexes)
        species_type = self.model.species_types.create(id='complex_70S_IA', name='complex_70S_IA', type=5)
        species = species_type.species.create(compartment=compartment)
        species.concentration = wc_lang.core.Concentration(value=3000, units='molecules')

        species_type = self.model.species_types.create(id='complex_70S_A', name='complex_70S_A', type=5)
        species = species_type.species.create(compartment=compartment)
        species.concentration = wc_lang.core.Concentration(value=3000, units='molecules')

        # Create both functional and afunctional forms (_att: attached to RNA) of every protein in KB
        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):

            # Add functional form of protein
            species_type = self.model.species_types.create(
                id                = protein.id,
                type              = wc_lang.SpeciesTypeType.protein,
                name              = protein.name,
                structure         = protein.get_seq(cds=False),
                empirical_formula = protein.get_empirical_formula(cds=False),
                molecular_weight  = protein.get_mol_wt(cds=False),
                charge            = protein.get_charge(cds=False))

            species = species_type.species.create(compartment=compartment)
            species.concentration = wc_lang.core.Concentration(value=100, units='molecules')

            # Add inactive form of protein, attached to ribosome
            species_type = self.model.species_types.create(
                id                = protein.id+'_att',
                type              = wc_lang.SpeciesTypeType.pseudo_species,
                name              = protein.name+'_att',
                structure         = protein.get_seq(cds=False),
                empirical_formula = protein.get_empirical_formula(cds=False),
                molecular_weight  = protein.get_mol_wt(cds=False),
                charge            = protein.get_charge(cds=False))

            species = species_type.species.create(compartment=compartment)
            species.concentration = wc_lang.core.Concentration(value=100, units=wc_lang.ConcentrationUnit.molecules)

    def gen_reactions(self):
        """ Generate a set of 3 reqactions (initation, elongation, termination) for each protein """

        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        amino_acids = {'S': ['tRNA-Ser', 'Rna_MPNt09', 'Rna_MPNt24', 'Rna_MPNt25', 'Rna_MPNt26', 'Rna_MPNt03'],
                       'L': ['tRNA-Leu', 'Rna_MPNt19', 'Rna_MPNt27', 'Rna_MPNt35', 'Rna_MPNt36'],
                       'R': ['tRNA-Arg', 'Rna_MPNt15', 'Rna_MPNt17', 'Rna_MPNt37'],
                       'T': ['tRNA-Thr', 'Rna_MPNt04',  'Rna_MPNt29', 'Rna_MPNt31'],
                       'G': ['tRNA-Gly', 'Rna_MPNt18', 'Rna_MPNt14'],
                       'F': ['tRNA-Phe', 'Rna_MPNt12',  'Rna_MPNt13'],
                       'W': ['tRNA-Trp', 'Rna_MPNt16',  'Rna_MPNt23'],
                       'K': ['tRNA-Lys', 'Rna_MPNt20', 'Rna_MPNt28'],
                       'I': ['tRNA-Ile', 'Rna_MPNt02',  'Rna_MPNt08'],
                       'A': ['tRNA-Ala', 'Rna_MPNt01'],
                       'M': ['tRNA-Met', 'Rna_MPNt07'],
                       'Q': ['tRNA-Gln', 'Rna_MPNt21'],
                       'E': ['tRNA-Glu', 'Rna_MPNt32'],
                       'P': ['tRNA-Pro', 'Rna_MPNt06'],
                       'V': ['tRNA-Val', 'Rna_MPNt30'],
                       'C': ['tRNA-Cys', 'Rna_MPNt05'],
                       'Y': ['tRNA-Tyr', 'Rna_MPNt22'],
                       'H': ['tRNA-His', 'Rna_MPNt34'],
                       'N': ['tRNA-Asn', 'Rna_MPNt33'],
                       'D': ['tRNA-Asp', 'Rna_MPNt11']}

        # Add translation initating reactions
        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):
            reaction = wc_lang.core.Reaction(id='translation_init_' + protein.id, submodel=submodel)
            reaction.participants=[]

            # Adding reaction participants LHS
            specie = self.model.species_types.get_one(id='complex_70S_IA').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))
            specie = self.model.species_types.get_one(id='GTP').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))
            specie = self.model.species_types.get_one(id='prot_MPN187').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))
            specie = self.model.species_types.get_one(id='prot_MPN631').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))
            specie = self.model.species_types.get_one(id='prot_MPN227').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))

            #specie = self.model.species_types.get_one(id='prot_MPN227').species.get_one(compartment=compartment)
            #reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1)

            # Adding reaction participants RHS
            specie = self.model.species_types.get_one(id='complex_70S_A').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))
            specie = self.model.species_types.get_one(id='GDP').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))
            specie = self.model.species_types.get_one(id='P').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))
            specie = self.model.species_types.get_one(id='prot_MPN187').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))
            specie = self.model.species_types.get_one(id='prot_MPN631').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))
            specie = self.model.species_types.get_one(id='prot_MPN227').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))

        # Add translation elongation reactions
        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):
            reaction = wc_lang.core.Reaction(id='translation_elon_' + protein.id, submodel=submodel)
            n_steps = len(protein.get_seq(cds=False))
            reaction.participants=[]

            # Adding reaction participants LHS
            specie = self.model.species_types.get_one(id='complex_70S_A').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))
            specie = self.model.species_types.get_one(id='prot_MPN665').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-n_steps))
            specie = self.model.species_types.get_one(id='prot_MPN631').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-n_steps))
            specie = self.model.species_types.get_one(id='prot_MPN227').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-n_steps))
            specie = self.model.species_types.get_one(id='GTP').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-2*n_steps))

            # tRNAs - add more testing to this bit
            for letter in list(set(protein.get_seq(cds=False).tostring())):
                if letter == '*':
                    continue

                n = protein.get_seq(cds=False).tostring().count(letter)
                specie = self.model.species_types.get_or_create(
                    id=amino_acids[letter][1], type=4).species.get_or_create(compartment=compartment)
                reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-n))

            # Adding reaction participants RHS
            specie = self.model.species_types.get_one(id=protein.id+'_att').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))
            specie = self.model.species_types.get_one(id='GDP').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=2*n_steps))
            specie = self.model.species_types.get_one(id='P').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=2*n_steps))

        # Add translation termination reactions
        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):
            reaction = wc_lang.core.Reaction(id='translation_term_' + protein.id, submodel=submodel)
            n_steps = len(protein.get_seq(cds=False))
            reaction.participants=[]

            # Adding reaction participants LHS
            specie = self.model.species_types.get_one(id='prot_MPN361').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))
            specie = self.model.species_types.get_one(id='GTP').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))
            specie = self.model.species_types.get_one(id=protein.id+'_att').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=-1))

            # Adding reaction participants RHS
            specie = self.model.species_types.get_one(id=protein.id).species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))
            specie = self.model.species_types.get_one(id='GDP').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))
            specie = self.model.species_types.get_one(id='P').species.get_one(compartment=compartment)
            reaction.participants.add(specie.species_coefficients.get_or_create(coefficient=1))

    def gen_rate_laws(self):
        """ Generate the rate laws associate dwith reactions """
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        for reaction in submodel.reactions:
            exp= 'k_cat'
            mod = []
            rate_eq = None

            for participant in reaction.participants:
                if participant.coefficient > 0:
                    continue

                if participant.coefficient < 0:
                    exp = exp + ' * (' + participant.species.id() + '/ (k_m + ' + participant.species.id() + '))'
                    mod.append(self.model.species_types.get_one(
                        id=participant.species.species_type.id).species.get_one(compartment=compartment))

            for rxn in submodel.reactions:
              for rl in rxn.rate_laws:
                  if rl.equation.expression == exp:
                      rate_eq = rl.equation

            if rate_eq is None:
              rate_eq = wc_lang.core.RateLawEquation(expression=exp, modifiers=mod)

            rate_law = wc_lang.core.RateLaw(reaction=reaction,
                                          direction=wc_lang.core.RateLawDirection.forward,
                                          equation=rate_eq,
                                          k_cat=10,
                                          k_m=1)
