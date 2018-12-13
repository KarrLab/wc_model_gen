""" Generator for transcription submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen
import wc_lang
import wc_kb
import numpy
import math


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
        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rna_kbs:

            rna_model = model.species_types.get_one(id=rna_kb.id).species.get_one(compartment=cytosol)
            reaction = model.reactions.get_or_create(submodel=submodel, id=rna_kb.id.replace('rna_', 'transcription_'))
            reaction.name = rna_kb.id.replace('rna_', 'transcription_')
            reaction.participants = []
            seq = rna_kb.get_seq()

            # Adding participants to LHS
            reaction.participants.add(atp.species_coefficients.get_or_create(coefficient=-seq.count('A')))
            reaction.participants.add(ctp.species_coefficients.get_or_create(coefficient=-seq.count('C')))
            reaction.participants.add(gtp.species_coefficients.get_or_create(coefficient=-seq.count('G')))
            reaction.participants.add(utp.species_coefficients.get_or_create(coefficient=-seq.count('U')))
            reaction.participants.add(h.species_coefficients.get_or_create(coefficient=-(rna_kb.get_len() - 1)))

            # Adding participants to RHS
            reaction.participants.add(rna_model.species_coefficients.get_or_create(coefficient=1))
            reaction.participants.add(ppi.species_coefficients.get_or_create(coefficient=rna_kb.get_len()))
            reaction.participants.add(h2o.species_coefficients.get_or_create(coefficient=rna_kb.get_len() - 1))

            # Add RNA polymerease
            for rnap_kb in cell.observables.get_one(id='rna_polymerase_obs').species:
                rnap_species_type_model = model.species_types.get_one(id=rnap_kb.species.species_type.id)
                rnap_model = rnap_species_type_model.species.get_one(compartment=cytosol)

                reaction.participants.add(rnap_model.species_coefficients.get_or_create(coefficient=(-1)*rnap_kb.coefficient))
                reaction.participants.add(rnap_model.species_coefficients.get_or_create(coefficient=rnap_kb.coefficient))

    def gen_phenom_rates(self):
        """ Generate rate laws with exponential dynamics """
        submodel = self.model.submodels.get_one(id='transcription')
        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        cell_cycle_len = self.knowledge_base.cell.properties.get_one(id='cell_cycle_len').value

        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):
            self.gen_phenom_rate_law_eq(specie_type_kb=rna_kb,
                                        reaction=reaction,
                                        half_life=rna_kb.half_life,
                                        cell_cycle_len=cell_cycle_len)

    def gen_mechanistic_rates(self):
        """ Generate rate laws with calibrated dynamics """
        submodel = self.model.submodels.get_one(id='transcription')
        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        cell_cycle_len = self.knowledge_base.cell.properties.get_one(id='cell_cycle_len').value

        for rna_kb, reaction in zip(rnas_kb, self.submodel.reactions):
            self.gen_mechanistic_rate_law_eq(specie_type_kb=rna_kb,
                                             submodel=submodel,
                                             reaction=reaction,
                                             beta=1.,
                                             half_life=rna_kb.half_life,
                                             cell_cycle_len=cell_cycle_len)
