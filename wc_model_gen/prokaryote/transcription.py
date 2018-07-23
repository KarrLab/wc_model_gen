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
        compartment = self.model.compartments.get_one(id='c')
        for rna in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType):

            species_type = self.model.species_types.create(
                id=rna.id,
                name=rna.name,
                structure=rna.get_seq(),
                empirical_formula=rna.get_empirical_formula(),
                molecular_weight=rna.get_mol_wt(),
                charge=rna.get_charge(),
                type=wc_lang.SpeciesTypeType.rna)

            species = species_type.species.create(compartment=compartment)
            species.concentration = wc_lang.core.Concentration(value=0.00001, units=wc_lang.ConcentrationUnit.uM)

    def gen_reactions(self):
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        atp = self.model.species_types.get_one(id='ATP').species.get_one(compartment=compartment)
        ctp = self.model.species_types.get_one(id='CTP').species.get_one(compartment=compartment)
        gtp = self.model.species_types.get_one(id='GTP').species.get_one(compartment=compartment)
        utp = self.model.species_types.get_one(id='UTP').species.get_one(compartment=compartment)
        h2o = self.model.species_types.get_one(id='H2O').species.get_one(compartment=compartment)
        ppi = self.model.species_types.get_one(id='PPI').species.get_one(compartment=compartment)
        h   = self.model.species_types.get_one(id='H').species.get_one(compartment=compartment)

        # Generate a reaction for each Rna
        for rna in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.RnaSpeciesType):

            rna_specie = self.model.species_types.get_one(id=rna.id).species.get_one(compartment=compartment)
            reaction = wc_lang.core.Reaction(id='transcription_' + rna.id, submodel=submodel)
            reaction.participants=[]
            seq = rna.get_seq()

            # Adding reaction participants LHS
            reaction.participants.add(atp.species_coefficients.get_or_create(coefficient=-seq.count('A')))
            reaction.participants.add(ctp.species_coefficients.get_or_create(coefficient=-seq.count('C')))
            reaction.participants.add(gtp.species_coefficients.get_or_create(coefficient=-seq.count('G')))
            reaction.participants.add(utp.species_coefficients.get_or_create(coefficient=-seq.count('U')))

            # Adding reaction participants RHS
            reaction.participants.add(rna_specie.species_coefficients.get_or_create(coefficient=1))
            reaction.participants.add(h2o.species_coefficients.get_or_create(coefficient=(rna.get_len()-1)))
            reaction.participants.add(ppi.species_coefficients.get_or_create(coefficient=len(seq)))
            reaction.participants.add(h.species_coefficients.get_or_create(coefficient=(rna.get_len()-1)))

    def gen_rate_laws(self):
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        for reaction in submodel.reactions:
            exp= 'k_cat'
            mod = []
            rate_eq = None

            for participant in reaction.participants:
                if participant.coefficient > 0:
                    continue

                if participant.coefficient < 0:
                    exp = exp + ' * (' + participant.species.id() + '/ (k_m + ' + participant.species.id() + '))'
                    mod.append(self.model.species_types.get_one(
                        id=participant.species.species_type.id).species.get_one(compartment=compartment))

            for rxn in submodel.reactions:
                for rl in rxn.rate_laws:
                    if rl.equation.expression == exp:
                        rate_eq = rl.equation

            if rate_eq is None:
                rate_eq = wc_lang.core.RateLawEquation(expression=exp, modifiers=mod)

            rate_law = wc_lang.core.RateLaw(reaction=reaction,
                                            direction=wc_lang.core.RateLawDirection.forward,
                                            equation=rate_eq,
                                            k_cat=1,
                                            k_m=1)
