""" Generating wc_lang formatted models from knowledge base.
:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT

TODO: aminoacyl-tRnas are missign from observables
"""

import wc_kb
import wc_lang
import wc_model_gen
import numpy
import scipy

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
        pi =  model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)
        initiation_factors = model.observables.get_one(id='translation_init_factors_obs').expression.species[0]
        elongation_factors = model.observables.get_one(id='translation_elongation_factors_obs').expression.species[0]
        release_factors = model.observables.get_one(id='translation_release_factors_obs').expression.species[0]

        bases = "TCAG"
        codons = [a + b + c for a in bases for b in bases for c in bases]

        proteins_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb in proteins_kbs:

            protein_model = model.species_types.get_one(id=protein_kb.id).species.get_one(compartment=cytosol)
            n_steps = protein_kb.get_len()
            rxn = submodel.reactions.get_or_create(id=protein_kb.id.replace('prot_', 'translation_'))
            rxn.participants = []

            # Adding participants to LHS
            rxn.participants.add(gtp.species_coefficients.get_or_create(coefficient=-(n_steps+2)))
            rxn.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=-n_steps))
            rxn.participants.add(release_factors.species_coefficients.get_or_create(coefficient=-1))

            # Add tRNAs to LHS
            for codon in codons:
                if codon not in ['TAG', 'TAA', 'TGA']:
                    n=0
                    for base in range(0,len(protein_kb.gene.get_seq()),3):
                        n += str(protein_kb.gene.get_seq()[base:base+3]).count(codon)
                    if n > 0:
                        trna = model.observables.get_one(id='tRNA_'+codon+'_obs').expression.species[0]
                        rxn.participants.add(trna.species_coefficients.get_or_create(coefficient=-n))

            # Adding participants to RHS
            # Taking care of special cases of 'self-production': there should only be one SCoeff not 2!
            if protein_model==initiation_factors:
                rxn.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=2))
                rxn.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps))
                rxn.participants.add(release_factors.species_coefficients.get_or_create(coefficient=1))

            elif protein_model==elongation_factors:
                rxn.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps+1))
                rxn.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=1))
                rxn.participants.add(release_factors.species_coefficients.get_or_create(coefficient=1))

            elif protein_model==release_factors:
                rxn.participants.add(release_factors.species_coefficients.get_or_create(coefficient=2))
                rxn.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=1))
                rxn.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps))

            else:
                rxn.participants.add(protein_model.species_coefficients.get_or_create(coefficient=1))
                rxn.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=1))
                rxn.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps))
                rxn.participants.add(release_factors.species_coefficients.get_or_create(coefficient=1))

            rxn.participants.add(gdp.species_coefficients.get_or_create(coefficient=n_steps+2))
            rxn.participants.add(pi.species_coefficients.get_or_create(coefficient=2*n_steps))

            # Add ribosome
            for ribosome_kb in cell.observables.get_one(id='ribosome_obs').species:
                ribosome_species_type_model = model.species_types.get_one(id=ribosome_kb.species.species_type.id)
                ribosome_model = ribosome_species_type_model.species.get_one(compartment=cytosol)

                rxn.participants.add(ribosome_model.species_coefficients.get_or_create(coefficient=(-1)*ribosome_kb.coefficient))
                rxn.participants.add(ribosome_model.species_coefficients.get_or_create(coefficient=ribosome_kb.coefficient))

    def gen_rate_laws(self):
        """ Choose dynamics for the model """

        rate_law_dynamics = self.options.get('rate_law_dynamics')
        if rate_law_dynamics=='exponential':
            self.gen_rate_laws_exp()

        elif rate_law_dynamics=='calibrated':
            self.gen_rate_laws_cal()

    def gen_phenomenological_rates(self):
        """ Generate rate laws with exponential dynamics """

        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='translation')
        cell_cycle_length = cell.properties.get_one(id='cell_cycle_length').value

        proteins_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb, rxn in zip(proteins_kbs, submodel.reactions):
            protein_model = model.species_types.get_one(id=protein_kb.id).species[0]

            if protein_kb.half_life == 0:
                protein_kb.half_life = 12*60*60

            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = '({} / {} + {} / {}) * {}'.format(numpy.log(2), protein_kb.half_life,
                                                           numpy.log(2), cell_cycle_length,
                                                           protein_model.id())

            rate_law.equation = wc_lang.RateLawEquation(expression = expression)
            rate_law.equation.modifiers.append(protein_model)

    def gen_mechanistic_rates(self):
        """ Generate rate laws associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='translation')
        mean_volume = cell.properties.get_one(id='initial_volume').value
        mean_cell_cycle_length = cell.properties.get_one(id='cell_cycle_length').value
        cytosol = cell.compartments.get_one(id='c')

        proteins_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb, rxn in zip(proteins_kbs, submodel.reactions):

            protein_model = model.species_types.get_one(id=protein_kb.id).species[0]
            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = 'k_cat*'
            modifiers = []
            rate_avg = ''
            beta = 2

            #TODO: replace with calculation of avg half life; 553s is avg of Mycoplasma RNAs
            if protein_kb.half_life == 0:
                protein_kb.half_life = 12*60*60

            for participant in rxn.participants:
                if participant.coefficient < 0:
                    avg_conc = (3/2)*participant.species.concentration.value
                    modifiers.append(participant.species)
                    rate_avg   += '({}/({}+({}*{})))*'.format(avg_conc, avg_conc, beta, avg_conc)
                    expression += '({}/({}+{}*{}))*'.format(participant.species.id(),
                                                              participant.species.id(),
                                                              beta,
                                                              participant.species.concentration.value)

            # Clip off trailing * character
            expression = expression[:-1]
            rate_avg = rate_avg[:-1]

            # Create / add rate law equation
            if 'rate_law_equation' not in locals():
                rate_law_equation = wc_lang.RateLawEquation(expression=expression, modifiers=modifiers)

            rate_law.equation = rate_law_equation

            # Calculate k_cat
            exp_expression = '({}*(1/{}+1/{})*{})'.format(
                                numpy.log(2),
                                cell.properties.get_one(id='cell_cycle_length').value,
                                protein_kb.half_life,
                                3/2*protein_kb.species.get_one(compartment=cytosol).concentrations.value)
                                #This should have units of M

            rate_law.k_cat = eval(exp_expression) / eval(rate_avg)
