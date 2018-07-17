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


class RnaDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for RNA degradation submodel """

    def gen_compartments(self):
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get_or_create(id='c')
        cytosol.name = 'cytosol'
        cytosol.initial_volume = cell.properties.get_one(
            id='mean_volume').value

    def gen_species(self):
        """ Generate species associated with submodel """
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

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
                rxn = submodel.reactions.get_or_create(id=kb_rna.id.replace('rna_', 'rna_degradation_'))
                rxn.name = kb_rna.name.replace('RNA ', 'RNA degradation ')
            else:
                rxn = submodel.reactions.get_or_create(id='rna_degradation_'+str(kb_rna.id))
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

        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna, rxn in zip(rnas, self.submodel.reactions):
            rl = rxn.rate_laws.create()
            rl.direction = wc_lang.RateLawDirection.forward
            rl.equation = wc_lang.RateLawEquation(
                expression='k_cat * {0}[c] / (k_m + {0}[c])'.format(rna.id))
            rl.k_cat = 2 * numpy.log(2) / rna.half_life
            rl.k_m = rna.concentration
            rl.equation.modifiers.append(rxn.participants[0].species)
