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

        # print('here')
        submodel = self.submodel
        cytosol = self.model.compartments.get_one(id='c')
        cell = self.knowledge_base.cell
        model = self.model

        # Initiate ribosome species types (complexes)
        species_type = self.model.species_types.create(
            id='complex_70S_IA', name='complex_70S_IA', type=wc_lang.SpeciesTypeType.pseudo_species)
        species_type.molecular_weight = 1  # placeholder
        species = species_type.species.create(compartment=cytosol)
        species.concentration = wc_lang.core.Concentration(
            value=1e-2, units=wc_lang.ConcentrationUnit.M)

        species_type = self.model.species_types.create(
            id='complex_70S_A', name='complex_70S_A', type=wc_lang.SpeciesTypeType.pseudo_species)
        species_type.molecular_weight = 1  # placeholder
        species = species_type.species.create(compartment=cytosol)
        species.concentration = wc_lang.core.Concentration(
            value=1e-2, units=wc_lang.ConcentrationUnit.M)

        # get or create RNA species
        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna in rnas:
            species_type = model.species_types.get_or_create(id=rna.id)
            if not species_type.name:
                species_type.name = rna.name
                species_type.type = wc_lang.SpeciesTypeType.rna
                species_type.structure = rna.get_seq()
                species_type.empirical_formula = rna.get_empirical_formula()
                species_type.molecular_weight = rna.get_mol_wt()
                species_type.charge = rna.get_charge()
                species = species_type.species.get_or_create(
                    compartment=cytosol)
                species.concentration = wc_lang.Concentration(
                    value=rna.concentration, units=wc_lang.ConcentrationUnit.M)

        # Create both functional and afunctional form (_att: attached to RNA) of every protein in KB
        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):

            species_type = self.model.species_types.get_or_create(id=protein.id)
            if not species_type.name:
                # Add functional form of protein
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


            # Add inactive form of protein, attached to ribosome
            species_type = self.model.species_types.create(
                id=protein.id+'_att',
                type=wc_lang.SpeciesTypeType.protein,
                name=protein.name+'_att',
                structure=protein.get_seq(),
                empirical_formula=protein.get_empirical_formula(),
                molecular_weight=protein.get_mol_wt(),
                charge=protein.get_charge())

            species = species_type.species.create(compartment=cytosol)
            species.concentration = wc_lang.core.Concentration(
                value=0, units=wc_lang.ConcentrationUnit.M)

        for kb_observable in self.knowledge_base.cell.observables:
            model_observable = self.model.observables.get_or_create(
                id=kb_observable.id)
            if not model_observable.name:
                model_observable.name = kb_observable.name
                for kb_species_coefficient in kb_observable.species:
                    kb_species = kb_species_coefficient.species
                    kb_species_type = kb_species.species_type
                    kb_compartment = kb_species.compartment
                    model_species_type = model.species_types.get_one(
                        id=kb_species_type.id)
                    model_species = model_species_type.species.get_one(
                        compartment=model.compartments.get_one(id=kb_compartment.id))
                    model_coefficient = kb_species_coefficient.coefficient
                    model_species_coefficient = wc_lang.SpeciesCoefficient()
                    model_species_coefficient.species = model_species
                    model_species_coefficient.coefficient = model_coefficient

                    model_observable.species.append(model_species_coefficient)

                for kb_observable_observable in kb_observable.observables:
                    model_observable_observable = model.observables.get_or_create(
                        id=kb_observable_observable.id)
                    model_observable.observables.append(
                        model_observable_observable)

    def gen_reactions(self):
        """ Generate a set of 3 reqactions (initation, elongation, termination) for each protein """
        # print('here')
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        prots = self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType)

        amino_acids = {}
        amino_acids['S'] = 'tRNA_Ser'
        amino_acids['L'] = 'tRNA_Leu'
        amino_acids['R'] = 'tRNA_Arg'
        amino_acids['T'] = 'tRNA_Thr'
        amino_acids['G'] = 'tRNA_Gly'
        amino_acids['F'] = 'tRNA_Phe'
        amino_acids['W'] = 'tRNA_Trp'
        amino_acids['K'] = 'tRNA_Lys'
        amino_acids['I'] = 'tRNA_Ile'
        amino_acids['A'] = 'tRNA_Ala'
        amino_acids['M'] = 'tRNA_Met'
        amino_acids['Q'] = 'tRNA_Gln'
        amino_acids['E'] = 'tRNA_Glu'
        amino_acids['P'] = 'tRNA_Pro'
        amino_acids['V'] = 'tRNA_Val'
        amino_acids['C'] = 'tRNA_Cys'
        amino_acids['Y'] = 'tRNA_Tyr'
        amino_acids['H'] = 'tRNA_His'
        amino_acids['N'] = 'tRNA_Asn'
        amino_acids['D'] = 'tRNA_Asp'

        # Add translation initating reactions
        for protein in prots:
                reaction = wc_lang.core.Reaction(
                    id='translation_init_' + protein.id, submodel=submodel)
                reaction.name = protein.id

                # Adding reaction participants LHS
                specie = self.model.species_types.get_one(
                    id='complex_70S_IA').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-1))
                specie = self.model.species_types.get_one(
                    id='gtp').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-1))

                # Adding reaction participants RHS
                specie = self.model.species_types.get_one(
                    id='complex_70S_A').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))
                specie = self.model.species_types.get_one(
                    id='gdp').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))
                specie = self.model.species_types.get_one(
                    id='pi').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))

        # Add translation elongation reactions
        for protein in prots:
                reaction = wc_lang.core.Reaction(
                    id='translation_elon_' + protein.id, submodel=submodel)
                reaction.name = protein.id
                n_steps = len(protein.get_seq())

                # Adding reaction participants LHS
                specie = self.model.species_types.get_one(
                    id='complex_70S_A').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-1))
                specie = self.model.species_types.get_one(
                    id='gtp').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-2 * n_steps))

                for letter in list(set(str(protein.get_seq()))):
                    if letter == '*':
                        continue

                    n = protein.get_seq().tostring().count(letter)
                    specie = self.model.observables.get_one(
                        id=amino_acids[letter]+'_obs').species[0].species
                    reaction.participants.add(
                        specie.species_coefficients.get_or_create(coefficient=-n))

                # Adding reaction participants RHS
                specie = self.model.species_types.get_one(
                    id=protein.id+'_att').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))
                specie = self.model.species_types.get_one(
                    id='gdp').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=2 * n_steps))
                specie = self.model.species_types.get_one(
                    id='pi').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=2 * n_steps))
                specie = self.model.species_types.get_one(
                    id='complex_70S_A').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))

        # Add translation termination reactions
        for protein in prots:
                reaction = wc_lang.core.Reaction(
                    id='translation_term_' + protein.id, submodel=submodel)
                reaction.name = protein.id
                n_steps = len(protein.get_seq())

                # Adding reaction participants LHS

                specie = self.model.species_types.get_one(
                    id='gtp').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-1))
                specie = self.model.species_types.get_one(
                    id=protein.id+'_att').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-1))
                specie = self.model.species_types.get_one(
                    id='complex_70S_A').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-1))

                # Adding reaction participants RHS
                specie = self.model.species_types.get_one(
                    id=protein.id).species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))
                specie = self.model.species_types.get_one(
                    id='gdp').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))
                specie = self.model.species_types.get_one(
                    id='pi').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))
                specie = self.model.species_types.get_one(
                    id='complex_70S_IA').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=1))

    def gen_rate_laws(self):
        """ Generate the rate laws associate dwith reactions """
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        IF = self.model.observables.get_one(
            id='IF_obs')
        IF = IF.species[0].species.species_type
        EF = self.model.observables.get_one(
            id='EF_obs')
        EF = EF.species[0].species.species_type
        RF = self.model.observables.get_one(
            id='RF_obs')
        RF = RF.species[0].species.species_type
        

        for reaction in submodel.reactions:
            exp = 'k_cat'
            mod = []

            rate_eq = None

            for participant in reaction.participants:
                if participant.coefficient > 0:
                    continue

                if participant.coefficient < 0:
                    exp = exp + ' * (' + participant.species.id() + \
                        '/ (k_m + ' + participant.species.id() + '))'
                    mod.append(self.model.species_types.get_one(
                        id=participant.species.species_type.id).species.get_one(compartment=compartment))

            if reaction.id.startswith('translation_init_'):
                        exp = exp + ' * ({}[c]'.format(IF.id) + \
                            '/ (k_m +{}[c]))'.format(IF.id)
                        mod.append(IF.species.get_one(compartment = compartment))
                    
            elif reaction.id.startswith('translation_elon_'):
                    exp = exp + ' * ({}[c]'.format(EF.id) + \
                        '/ (k_m +{}[c]))'.format(EF.id)
                    mod.append(EF.species.get_one(compartment = compartment))

            else:
                    exp = exp + ' * ({}[c]'.format(RF.id) + \
                        '/ (k_m +{}[c]))'.format(RF.id)
                    mod.append(RF.species.get_one(compartment = compartment))
                

            for rxn in submodel.reactions:
                for rl in rxn.rate_laws:
                    if rl.equation.expression == exp:
                        rate_eq = rl.equation

            if rate_eq is None:
                rate_eq = wc_lang.core.RateLawEquation(
                    expression=exp, modifiers=mod)

            rate_law = wc_lang.core.RateLaw(reaction=reaction,
                                            direction=wc_lang.core.RateLawDirection.forward,
                                            equation=rate_eq,
                                            k_cat=1,
                                            k_m=1)
