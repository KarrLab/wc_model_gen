""" Generator for transcription submodels based on KBs for random in silico organisms

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

class TranscriptionSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for transcription submodel """

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        submodel = self.submodel
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        # Get species involved in reaction
        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        # Create reaction for each RNA
        # TODO: could eliminate refering to rna_kb: rna_model already has all species
        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rna_kbs:

            rna_model = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            seq = rna_kb.get_seq()
            rxn = submodel.reactions.get_or_create(id=rna_kb.id.replace('rna_', 'transcription_'))
            rxn.participants = []

            # Adding participants to LHS
            rxn.participants.add(atp.species_coefficients.get_or_create(coefficient=-seq.count('A')))
            rxn.participants.add(ctp.species_coefficients.get_or_create(coefficient=-seq.count('C')))
            rxn.participants.add(gtp.species_coefficients.get_or_create(coefficient=-seq.count('G')))
            rxn.participants.add(utp.species_coefficients.get_or_create(coefficient=-seq.count('U')))
            rxn.participants.add(h.species_coefficients.get_or_create(coefficient=-(rna_kb.get_len() - 1)))

            # Adding participants to RHS
            rxn.participants.add(rna_model.species_coefficients.get_or_create(coefficient=1))
            rxn.participants.add(ppi.species_coefficients.get_or_create(coefficient=rna_kb.get_len()))
            rxn.participants.add(h2o.species_coefficients.get_or_create(coefficient=rna_kb.get_len() - 1))

            # Add RNA polymerease
            for rnap_kb in cell.observables.get_one(id='rna_polymerase_obs').species:
                rnap_species_type_model = model.species_types.get_one(id=rnap_kb.species.species_type.id)
                rnap_model = rnap_species_type_model.species.get_one(compartment=cytosol)

                rxn.participants.add(rnap_model.species_coefficients.get_or_create(coefficient=(-1)*rnap_kb.coefficient))
                rxn.participants.add(rnap_model.species_coefficients.get_or_create(coefficient=rnap_kb.coefficient))

    def gen_phenomenological_rates(self):
        """ Generate rate laws with exponential dynamics """

        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='transcription')
        cell_cycle_length = cell.properties.get_one(id='cell_cycle_length').value

        rnas = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb, rxn in zip(rnas, self.submodel.reactions):

            rna_model = model.species_types.get_one(id=rna_kb.id).species[0]

            if rna_kb.half_life == 0:
                rna_kb.half_life = 553

            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = '({} / {} + {} / {}) * {}'.format(numpy.log(2), rna_kb.half_life,
                                                            numpy.log(2), cell_cycle_length,
                                                            rna_model.id())

            rate_law.equation = wc_lang.RateLawEquation(expression = expression)
            rate_law.equation.modifiers.append(rna_model)

    def gen_mechanistic_rates(self):
        """ Generate rate laws with calibrated dynamics """

        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='transcription')
        mean_volume = cell.properties.get_one(id='initial_volume').value
        mean_cell_cycle_length = cell.properties.get_one(id='cell_cycle_length').value

        rnas = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
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
                    avg_conc = participant.species.concentration.value# *(3/2)
                    modifiers.append(participant.species)
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
