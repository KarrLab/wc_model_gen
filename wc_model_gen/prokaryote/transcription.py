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
        cytosol = model.compartments.get_one(id='c')
        atp = model.species_types.get_one(
            id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(
            id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(
            id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(
            id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(
            id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(
            id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(
            id='h').species.get_one(compartment=cytosol)

        cell = self.knowledge_base.cell
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

            model_rna = model.species_types.get_one(
                id=kb_rna.id).species.get_one(compartment=cytosol)
            seq = kb_rna.get_seq()
            rxn.participants = []
            rxn.participants.add(atp.species_coefficients.get_or_create(
                coefficient=-seq.count('A')))
            rxn.participants.add(ctp.species_coefficients.get_or_create(
                coefficient=-seq.count('C')))
            rxn.participants.add(gtp.species_coefficients.get_or_create(
                coefficient=-seq.count('G')))
            rxn.participants.add(utp.species_coefficients.get_or_create(
                coefficient=-seq.count('U')))
            rxn.participants.add(h.species_coefficients.get_or_create(
                coefficient=-(kb_rna.get_len() - 1)))
            rxn.participants.add(
                model_rna.species_coefficients.get_or_create(coefficient=1))
            rxn.participants.add(ppi.species_coefficients.get_or_create(
                coefficient=kb_rna.get_len()))
            rxn.participants.add(h2o.species_coefficients.get_or_create(
                coefficient=kb_rna.get_len() - 1))

    def gen_rate_laws(self):
        """ Generate rate laws associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='transcription')

        mean_volume = cell.properties.get_one(id='initial_volume').value
        mean_doubling_time = cell.properties.get_one(id='doubling_time').value

        # http://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=2&id=106199
        poly_avg_conc = 3000/scipy.constants.Avogadro / cytosol.initial_volume
        rna_polymerase = model.observables.get_one(id='rna_polymerase_obs')

        for reaction in submodel.reactions:
            rate_law = reaction.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            rate_law_equation = wc_lang.RateLawEquation()

            rate_law_equation.modifiers.append(rna_polymerase.expression.species[0])
            expression = 'k_cat*({}/(k_m + {})'.format(rna_polymerase.expression.species[0].get_id(),
                                                       rna_polymerase.expression.species[0].get_id())

            for participant in reaction.participants:
                if participant.coefficient < 0:
                    rate_law_equation.modifiers.append(participant.species)
                    expression += '*({}/(k_m+{})'.format(
                                    participant.species.id(),
                                    participant.species.id())

            rate_law_equation.expression = expression
            rate_law.equation = rate_law_equation
