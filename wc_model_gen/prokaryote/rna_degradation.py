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
from wc_model_gen.prokaryote.species import SpeciesGenerator


class RnaDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for RNA degradation submodel """

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
        amp = model.species_types.get_one(
            id='amp').species.get_one(compartment=cytosol)
        cmp = model.species_types.get_one(
            id='cmp').species.get_one(compartment=cytosol)
        gmp = model.species_types.get_one(
            id='gmp').species.get_one(compartment=cytosol)
        ump = model.species_types.get_one(
            id='ump').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(
            id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(
            id='h').species.get_one(compartment=cytosol)

        cell = self.knowledge_base.cell
        kb_rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for kb_rna in kb_rnas:
            if kb_rna.id.startswith('rna_'):
                rxn = submodel.reactions.get_or_create(
                    id=kb_rna.id.replace('rna_', 'rna_degradation_'))
                rxn.name = kb_rna.name.replace('RNA ', 'RNA degradation ')
            else:
                rxn = submodel.reactions.get_or_create(
                    id='rna_degradation_'+str(kb_rna.id))
                rxn.name = 'RNA degradation '+str(kb_rna.name)

            model_rna = model.species_types.get_one(
                id=kb_rna.id).species.get_one(compartment=cytosol)
            seq = kb_rna.get_seq()
            rxn.participants = []
            rxn.participants.add(
                model_rna.species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(h2o.species_coefficients.get_or_create(
                coefficient=-(kb_rna.get_len() - 1)))
            rxn.participants.add(amp.species_coefficients.get_or_create(
                coefficient=seq.count('A')))
            rxn.participants.add(cmp.species_coefficients.get_or_create(
                coefficient=seq.count('C')))
            rxn.participants.add(gmp.species_coefficients.get_or_create(
                coefficient=seq.count('G')))
            rxn.participants.add(ump.species_coefficients.get_or_create(
                coefficient=seq.count('U')))
            rxn.participants.add(h.species_coefficients.get_or_create(
                coefficient=kb_rna.get_len() - 1))

    def gen_rate_laws(self):
        """ Generate rate laws associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        # http://bionumbers.hms.harvard.edu/bionumber.aspx?id=108959&ver=1&trm=average%20rnase%20concentration&org=
        deg_avg_conc = 5000/scipy.constants.Avogadro / cytosol.initial_volume

        deg_rnase = model.observables.get_one(id='deg_rnase_obs')

        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna, rxn in zip(rnas, self.submodel.reactions):
            rl = rxn.rate_laws.create()
            rl.direction = wc_lang.RateLawDirection.forward
            rl.equation = wc_lang.RateLawEquation(
                expression='{0}[c] * (((k_cat * {1}) / (k_m + {1})) + {2})'.format(rna.id, deg_rnase.species[0].species.id(), '0.1'))
            rl.k_cat = 2 * numpy.log(2) / rna.half_life
            rl.k_m = deg_avg_conc
            rl.equation.modifiers.append(deg_rnase.species[0].species)
            rl.equation.modifiers.append(rxn.participants[0].species)
