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
from wc_model_gen.prokaryote.species import SpeciesGenerator


class TranscriptionSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for transcription submodel """

    def gen_compartments(self):
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get_or_create(id='c')
        cytosol.name = 'cytosol'
        cytosol.initial_volume = cell.properties.get_one(
            id='initial_volume').value

    def gen_species(self):
        """ Generate species associated with submodel """

        speciesGen = SpeciesGenerator(self.knowledge_base, self.model)
        speciesGen.run()

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        submodel = self.submodel
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        rna_polymerase = model.observables.get_one(id='rna_polymerase_obs').expression.species[0]
        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        kb_rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for kb_rna in kb_rnas:
            if kb_rna.id.startswith('rna_'):
                rxn = submodel.reactions.get_or_create(
                    id=kb_rna.id.replace('rna_', 'transcription_'))
                rxn.name = kb_rna.name.replace('RNA ', 'Transcription')
            else:
                rxn = submodel.reactions.get_or_create(
                    id='transcription_'+str(kb_rna.id))
                rxn.name = 'Transcription '+str(kb_rna.name)

            model_rna = model.species_types.get_one(id=kb_rna.id).species.get_one(compartment=cytosol)
            seq = kb_rna.get_seq()
            rxn.participants = []

            rxn.participants.add(rna_polymerase.species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(atp.species_coefficients.get_or_create(coefficient=-seq.count('A')))
            rxn.participants.add(ctp.species_coefficients.get_or_create(coefficient=-seq.count('C')))
            rxn.participants.add(gtp.species_coefficients.get_or_create(coefficient=-seq.count('G')))
            rxn.participants.add(utp.species_coefficients.get_or_create(coefficient=-seq.count('U')))
            rxn.participants.add(h.species_coefficients.get_or_create(coefficient=-(kb_rna.get_len() - 1)))

            rxn.participants.add(rna_polymerase.species_coefficients.get_or_create(coefficient=1))
            rxn.participants.add(model_rna.species_coefficients.get_or_create(coefficient=1))
            rxn.participants.add(ppi.species_coefficients.get_or_create(coefficient=kb_rna.get_len()))
            rxn.participants.add(h2o.species_coefficients.get_or_create(coefficient=kb_rna.get_len() - 1))

    def gen_rate_laws(self):
        """ Generate rate laws for exponential submodel
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='transcription')
        cell_cycle_length = model.parameters.get_one(id='cell_cycle_length').value

        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna_kb, rxn in zip(rnas, self.submodel.reactions):

            rna_model = model.species_types.get_one(id=rna_kb.id).species[0]
            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = '({} / {} + {} / {}) * {}'.format(numpy.log(2), rna_kb.half_life,
                                                            numpy.log(2), cell_cycle_length,
                                                            rna_model.id())

            rate_law.equation = wc_lang.RateLawEquation(expression = expression)
            rate_law.equation.modifiers.append(rna_model)
        """

        """ Generate calibrated model """
        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='transcription')
        mean_volume = cell.properties.get_one(id='initial_volume').value
        mean_doubling_time = cell.properties.get_one(id='doubling_time').value
        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)

        for rna_kb, reaction in zip(rnas, submodel.reactions):

            rna_model = model.species_types.get_one(id=rna_kb.id).species[0]
            rate_law = reaction.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = 'k_cat*'
            modifiers = []
            rate_avg = ''
            beta = 2

            if rna_kb.half_life == 0:
                #TODO: replace with calculation of avg half life; 553s is avg of Mycoplasma RNAs
                rna_kb.half_life = 553

            for participant in reaction.participants:
                if participant.coefficient < 0:
                    avg_conc = (3/2)*participant.species.concentration.value
                    rate_avg += '({}/({}+({}*{})))*'.format(avg_conc, avg_conc, beta, avg_conc)
                    modifiers.append(participant.species)
                    expression += '({}/({}+(3/2)*{}))*'.format(participant.species.id(),
                                                              participant.species.id(),
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
                                model.parameters.get_one(id='cellCycleLength').value,
                                rna_kb.half_life,
                                3/2*rna_kb.concentration) #This should have units of M

            rate_law.k_cat = eval(exp_expression) / eval(rate_avg)
