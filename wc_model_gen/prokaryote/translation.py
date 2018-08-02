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
from wc_model_gen.prokaryote.species import SpeciesGenerator


class TranslationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generate translation submodel. """

    def gen_species(self):
        """ Generate a set of 2 protein species (_att: attached to ribosome and acitve form) for each protein """

        speciesGen = SpeciesGenerator(self.knowledge_base, self.model)
        speciesGen.run()
        
        # print('here')
        submodel = self.submodel
        cytosol = self.model.compartments.get_one(id='c')
        cell = self.knowledge_base.cell
        model = self.model

        # Create both functional and afunctional form (_att: attached to RNA) of every protein in KB
        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):

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
                specie = self.model.observables.get_one(
                    id='complex_70S_IA_obs').species[0].species
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-1))
                specie = self.model.species_types.get_one(
                    id='gtp').species.get_one(compartment=compartment)
                reaction.participants.add(
                    specie.species_coefficients.get_or_create(coefficient=-1))

                # Adding reaction participants RHS
                specie = self.model.observables.get_one(
                    id='complex_70S_A_obs').species[0].species
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


            
