""" Generating wc_lang formatted models from knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_kb
import wc_lang
import wc_model_gen


class TranscriptionSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generate transcription submodel. """

    def gen_species(self):
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        for rna in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType):
            species_type = self.model.species_types.create(
                id=rna.id,
                name=rna.name,
                structure=rna.get_seq(),
                empirical_formula=str(rna.get_empirical_formula()),
                molecular_weight=str(rna.get_mol_wt()),
                charge=str(rna.get_charge()),
                type=4)

            species = species_type.species.create(compartment=compartment)
            species.concentration = wc_lang.core.Concentration(value=5e-3, units='M')  # Default value

    def gen_reactions(self):
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        # Generate a reaction for each Rna
        for rna in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType):

            # Reaction attributes: 'id', 'submodel', 'rate_laws', 'participants', 'reversible', 'min_flux', 'max_flux'
            reaction = wc_lang.core.Reaction(id='transcription_' + rna.id, submodel=submodel)

            # Adding reaction participants LHS
            reaction.participants.create(species=self.model.species_types.get_one(id='ATP').species.get_one(compartment=compartment),
                                         coefficient=-(rna.get_seq().count('A')))
            reaction.participants.create(species=self.model.species_types.get_one(id='CTP').species.get_one(compartment=compartment),
                                         coefficient=-(rna.get_seq().count('C')))
            reaction.participants.create(species=self.model.species_types.get_one(id='GTP').species.get_one(compartment=compartment),
                                         coefficient=-(rna.get_seq().count('G')))
            reaction.participants.create(species=self.model.species_types.get_one(id='UTP').species.get_one(compartment=compartment),
                                         coefficient=-(rna.get_seq().count('U')))

            # Adding reaction participants RHS
            reaction.participants.create(species=self.model.species_types.get_one(id=rna.id).species.get_one(compartment=compartment),
                                         coefficient=1)
            reaction.participants.create(species=self.model.species_types.get_one(id='H2O').species.get_one(compartment=compartment),
                                         coefficient=1)
            reaction.participants.create(species=self.model.species_types.get_one(id='PPI').species.get_one(compartment=compartment),
                                         coefficient=len(rna.get_seq()))
            reaction.participants.create(species=self.model.species_types.get_one(id='H').species.get_one(compartment=compartment),
                                         coefficient=1)

    def gen_rate_laws(self):
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        for reaction in submodel.reactions:
            rate_eq = wc_lang.core.RateLawEquation(
                expression='k_cat'
                ' * (ATP[c] / (k_m + ATP[c]))'
                ' * (CTP[c] / (k_m + CTP[c]))'
                ' * (GTP[c] / (k_m + GTP[c]))'
                ' * (UTP[c] / (k_m + UTP[c]))',

                modifiers=[self.model.species_types.get_one(id='ATP').species.get_one(compartment=compartment),
                           self.model.species_types.get_one(id='CTP').species.get_one(compartment=compartment),
                           self.model.species_types.get_one(id='GTP').species.get_one(compartment=compartment),
                           self.model.species_types.get_one(id='UTP').species.get_one(compartment=compartment)])

            rate_law = wc_lang.core.RateLaw(reaction=reaction,
                                            direction=wc_lang.core.RateLawDirection.forward,
                                            equation=rate_eq,
                                            k_cat=1,
                                            k_m=1)
