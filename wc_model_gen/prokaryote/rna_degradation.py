""" Generator for RNA degradation submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import numpy
import scipy
import wc_kb
import wc_lang
import wc_model_gen
#from wc_model_gen.prokaryote.species import SpeciesGenerator

class RnaDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for RNA degradation submodel """

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='rna_degradation')
        cytosol = model.compartments.get_one(id='c')

        amp = model.species_types.get_one(id='amp').species.get_one(compartment=cytosol)
        cmp = model.species_types.get_one(id='cmp').species.get_one(compartment=cytosol)
        gmp = model.species_types.get_one(id='gmp').species.get_one(compartment=cytosol)
        ump = model.species_types.get_one(id='ump').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rna_kbs:

            rna_model = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            seq = rna_kb.get_seq()
            rxn = submodel.reactions.get_or_create(id=rna_kb.id.replace('rna_', 'degradation_'))
            rxn.participants = []

            # Adding participants to LHS
            rxn.participants.add(rna_model.species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(h2o.species_coefficients.get_or_create(coefficient=-(rna_kb.get_len() - 1)))

            # Adding participants to RHS
            rxn.participants.add(amp.species_coefficients.get_or_create(coefficient=seq.count('A')))
            rxn.participants.add(cmp.species_coefficients.get_or_create(coefficient=seq.count('C')))
            rxn.participants.add(gmp.species_coefficients.get_or_create(coefficient=seq.count('G')))
            rxn.participants.add(ump.species_coefficients.get_or_create(coefficient=seq.count('U')))
            rxn.participants.add(h.species_coefficients.get_or_create(coefficient=rna_kb.get_len() - 1))

            # Add members of the degradosome
            # Counterintuitively .specie is a KB species_coefficient object
            for degradosome_kb in cell.observables.get_one(id='degrade_rnase_obs').species:
                degradosome_species_type_model = model.species_types.get_one(id=degradosome_kb.species.species_type.id)
                degradosome_species_model = degradosome_species_type_model.species.get_one(compartment=cytosol)

                rxn.participants.add(degradosome_species_model.species_coefficients.get_or_create(coefficient=(-1)*degradosome_kb.coefficient))
                rxn.participants.add(degradosome_species_model.species_coefficients.get_or_create(coefficient=degradosome_kb.coefficient))

    def gen_phenomenological_rates(self):
        """ Generate rate laws with exponential dynamics """

        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='rna_degradation')

        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna_kb, rxn in zip(rnas, submodel.reactions):
            rna_model = model.species_types.get_one(id=rna_kb.id).species[0]
            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = '({} / {}) * {}'.format(numpy.log(2), rna_kb.half_life, rna_model.id())
            rate_law.equation = wc_lang.RateLawEquation(expression = expression)
            rate_law.equation.modifiers.append(rxn.participants[0].species)

    def gen_mechanistic_rates(self):
        """ Generate rate laws with calibrated dynamics """

        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='rna_degradation')
        mean_volume = cell.properties.get_one(id='initial_volume').value
        mean_cell_cycle_length = cell.properties.get_one(id='cell_cycle_length').value

        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna_kb, reaction in zip(rnas, submodel.reactions):

            rna_model = model.species_types.get_one(id=rna_kb.id).species[0]
            rate_law = reaction.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = 'k_cat*'
            modifiers = []
            rate_avg = ''
            beta = 1

            #TODO: replace with calculation of avg half life; 553s is avg of Mycoplasma RNAs
            if rna_kb.half_life == 0:
                rna_kb.half_life = 553

            for participant in reaction.participants:
                if participant.coefficient < 0:
                    modifiers.append(participant.species)
                    avg_conc = participant.species.concentration.value
                    rate_avg += '({}/({}+({}*{})))*'.format(avg_conc, avg_conc, beta, avg_conc)
                    expression += '({}/({}+({}*{})))*'.format(participant.species.id(),
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
                                rna_kb.half_life,
                                3/2*rna_kb.concentration) #This should have units of M

            rate_law.k_cat = eval(exp_expression) / eval(rate_avg)


        """
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        # http://bionumbers.hms.harvard.edu/bionumber.aspx?id=108959&ver=1&trm=average%20rnase%20concentration&org=
        deg_avg_conc = 5000/scipy.constants.Avogadro / cytosol.initial_volume

        deg_rnase = model.observables.get_one(
            id='degrade_rnase_obs').expression.species[0]

        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna, rxn in zip(rnas, self.submodel.reactions):
            rl = rxn.rate_laws.create()
            rl.direction = wc_lang.RateLawDirection.forward
            rl.equation = wc_lang.RateLawEquation(
                expression='{0}[c] * (((k_cat * {1}) / (k_m + {1})) + {2})'.format(rna.id, deg_rnase.id(), '0.1'))
            rl.k_cat = 2 * numpy.log(2) / rna.half_life
            rl.k_m = deg_avg_conc
            rl.equation.modifiers.append(deg_rnase)
            rl.equation.modifiers.append(rxn.participants[0].species)
        """
