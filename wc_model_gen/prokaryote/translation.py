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
import numpy


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
                    model_species_coefficient = model_species.species_coefficients.get_or_create(coefficient = model_coefficient)
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

        # Add translation initating reactions
        for protein in prots:
                reaction = submodel.reactions.get_or_create(id='translation_init_' + protein.id)
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
                reaction = submodel.reactions.get_or_create(id='translation_elon_' + protein.id)

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
                bases = "TCAG"
                codons = [a + b + c for a in bases for b in bases for c in bases]
                
                for codon in codons:
                    if codon not in ['TAG', 'TAA', 'TGA']:
                        n = str(protein.gene.get_seq()).count(codon)
                        if n > 0:
                            obs = self.model.observables.get_one(
                                id='tRNA_'+codon+'_obs')
                            specie = obs.species[0].species
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
                reaction = submodel.reactions.get_or_create(id='translation_term_' + protein.id)
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
        cell = self.knowledge_base.cell

        mean_volume = cell.properties.get_one(id='mean_volume').value
        mean_doubling_time = cell.properties.get_one(id='mean_doubling_time').value

        IF = self.model.observables.get_one(
            id='IF_obs')
        IF_avg_conc = 0.05 #placeholder
        EF = self.model.observables.get_one(
            id='EF_obs')
        EF_avg_conc = 0.05 #placeholder
        RF = self.model.observables.get_one(
            id='RF_obs')
        RF_avg_conc = 0.05 #placeholder

        exp = 'k_cat'
        
        init_eq = wc_lang.core.RateLawEquation(expression = exp + ' * ({}'.format(IF.id) + \
                  '/ (k_m +{}))'.format(IF.id))
        init_eq.observables.append(IF)

        elon_eq = wc_lang.core.RateLawEquation(expression = exp + ' * ({}'.format(EF.id) + \
                  '/ (k_m +{}))'.format(EF.id))
        elon_eq.observables.append(EF)


        term_eq = wc_lang.core.RateLawEquation(expression = exp + ' * ({}'.format(RF.id) + \
                                               '/ (k_m +{}))'.format(RF.id))
        term_eq.observables.append(RF)

        for reaction in self.submodel.reactions:
            if reaction.id.startswith('translation_init_'):
                        rl = reaction.rate_laws.create()
                        rl.direction = wc_lang.RateLawDirection.forward
                        prot_id = reaction.id[reaction.id.find('init_')+5:]
                        prot = cell.species_types.get_one(id=prot_id)
                        rl.k_cat = 2 * (numpy.log(2) / prot.half_life + numpy.log(2) / mean_doubling_time)
                        rl.equation = init_eq
                        rl.k_m = IF_avg_conc
            elif reaction.id.startswith('translation_elon_'):
                        rl = reaction.rate_laws.create()
                        rl.direction = wc_lang.RateLawDirection.forward
                        prot_id = reaction.id[reaction.id.find('elon_')+5:]
                        prot = cell.species_types.get_one(id=prot_id)
                        rl.k_cat = 2 * (numpy.log(2) / prot.half_life + numpy.log(2) / mean_doubling_time)
                        rl.equation = elon_eq
                        rl.k_m = EF_avg_conc

            else:
                        rl = reaction.rate_laws.create()
                        rl.direction = wc_lang.RateLawDirection.forward
                        prot_id = reaction.id[reaction.id.find('term_')+5:]
                        prot = cell.species_types.get_one(id=prot_id)
                        rl.k_cat = 2 * (numpy.log(2) / prot.half_life + numpy.log(2) / mean_doubling_time)
                        rl.equation = term_eq
                        rl.k_m = RF_avg_conc


            
