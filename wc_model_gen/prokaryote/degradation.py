""" Generating wc_lang formatted models from knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_lang
import wc_model_gen

class DegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):

    def gen_species(self):
        pass

    def gen_reactions(self):
        compartment = self.model.compartments.get_one(id='c')
        submodel = self.submodel
        model = self.model
        kb = self.knowledge_base

        # Get NMP molecules
        amp = model.species_types.get_one(id='AMP').species.get_one(compartment=compartment)
        cmp = model.species_types.get_one(id='CMP').species.get_one(compartment=compartment)
        gmp = model.species_types.get_one(id='GMP').species.get_one(compartment=compartment)
        ump = model.species_types.get_one(id='UMP').species.get_one(compartment=compartment)

        # Get AA molecules
        ala = model.species_types.get_one(id='Ala').species.get_one(compartment=compartment)
        arg = model.species_types.get_one(id='Arg').species.get_one(compartment=compartment)
        asn = model.species_types.get_one(id='Asn').species.get_one(compartment=compartment)
        asp = model.species_types.get_one(id='Asp').species.get_one(compartment=compartment)
        cys = model.species_types.get_one(id='Cys').species.get_one(compartment=compartment)
        gln = model.species_types.get_one(id='Gln').species.get_one(compartment=compartment)
        glu = model.species_types.get_one(id='Glu').species.get_one(compartment=compartment)
        gly = model.species_types.get_one(id='Gly').species.get_one(compartment=compartment)
        his = model.species_types.get_one(id='His').species.get_one(compartment=compartment)
        ile = model.species_types.get_one(id='Ile').species.get_one(compartment=compartment)
        leu = model.species_types.get_one(id='Leu').species.get_one(compartment=compartment)
        lys = model.species_types.get_one(id='Lys').species.get_one(compartment=compartment)
        met = model.species_types.get_one(id='Met').species.get_one(compartment=compartment)
        phe = model.species_types.get_one(id='Phe').species.get_one(compartment=compartment)
        pro = model.species_types.get_one(id='Pro').species.get_one(compartment=compartment)
        ser = model.species_types.get_one(id='Ser').species.get_one(compartment=compartment)
        thr = model.species_types.get_one(id='Thr').species.get_one(compartment=compartment)
        trp = model.species_types.get_one(id='Trp').species.get_one(compartment=compartment)
        tyr = model.species_types.get_one(id='Tyr').species.get_one(compartment=compartment)
        val = model.species_types.get_one(id='Val').species.get_one(compartment=compartment)

        # Loop through RNAs
        for rna_specie_type_model in model.species_types.get(type=4):
            # RNA species are represnted both in wc_lang and wc_kb, these are represented by
            # '_model' and '_kb' affixes respectively

            rna_specie_type_kb = kb.cell.species_types.get_one(id=rna_specie_type_model.id)
            rna_specie_model = model.species_types.get_one(id=rna_specie_type_model.id).species.get_one(compartment=compartment)

            reaction = wc_lang.core.Reaction(id='degradation_' + rna_specie_type_model.id, submodel=submodel)
            seq = rna_specie_type_kb.get_seq()
            reaction.participants=[]

            # Adding reaction participants LHS
            reaction.participants.add(rna_specie_model.species_coefficients.get_or_create(coefficient=-1))

            # Adding reaction participants RHS
            reaction.participants.add(amp.species_coefficients.get_or_create(coefficient=seq.count('A')))
            reaction.participants.add(cmp.species_coefficients.get_or_create(coefficient=seq.count('C')))
            reaction.participants.add(gmp.species_coefficients.get_or_create(coefficient=seq.count('G')))
            reaction.participants.add(ump.species_coefficients.get_or_create(coefficient=seq.count('U')))

        # Loop through proteins
        for protein_specie_type_model in self.model.species_types.get(type=2):

            protein_species_type_kb = kb.cell.species_types.get_one(id=protein_specie_type_model.id)
            #protein_specie_model = model.species_types.get_one(id=protein_specie_type_model.id).species.get_one(compartment=compartment)

            for protein_specie_model in protein_specie_type_model.species:
                seq = protein_species_type_kb.get_seq(cds=False)
                compartment = protein_specie_model.compartment

                reaction = wc_lang.core.Reaction(id='degradation_' +  protein_specie_model.species_type.id, submodel=submodel)
                reaction.participants=[]

                # Adding reaction participants LHS
                reaction.participants.add(protein_specie_model.species_coefficients.get_or_create(coefficient=-1))

                # Adding reaction participants RHS
                reaction.participants.add(ala.species_coefficients.get_or_create(coefficient=seq.count('A')))
                reaction.participants.add(arg.species_coefficients.get_or_create(coefficient=seq.count('R')))
                reaction.participants.add(asn.species_coefficients.get_or_create(coefficient=seq.count('N')))
                reaction.participants.add(asp.species_coefficients.get_or_create(coefficient=seq.count('D')))
                reaction.participants.add(cys.species_coefficients.get_or_create(coefficient=seq.count('C')))
                reaction.participants.add(gln.species_coefficients.get_or_create(coefficient=seq.count('Q')))
                reaction.participants.add(glu.species_coefficients.get_or_create(coefficient=seq.count('E')))
                reaction.participants.add(gly.species_coefficients.get_or_create(coefficient=seq.count('G')))
                reaction.participants.add(his.species_coefficients.get_or_create(coefficient=seq.count('H')))
                reaction.participants.add(ile.species_coefficients.get_or_create(coefficient=seq.count('I')))
                reaction.participants.add(leu.species_coefficients.get_or_create(coefficient=seq.count('L')))
                reaction.participants.add(lys.species_coefficients.get_or_create(coefficient=seq.count('K')))
                reaction.participants.add(met.species_coefficients.get_or_create(coefficient=seq.count('M')))
                reaction.participants.add(phe.species_coefficients.get_or_create(coefficient=seq.count('F')))
                reaction.participants.add(pro.species_coefficients.get_or_create(coefficient=seq.count('P')))
                reaction.participants.add(ser.species_coefficients.get_or_create(coefficient=seq.count('S')))
                reaction.participants.add(thr.species_coefficients.get_or_create(coefficient=seq.count('T')))
                reaction.participants.add(trp.species_coefficients.get_or_create(coefficient=seq.count('W')))
                reaction.participants.add(tyr.species_coefficients.get_or_create(coefficient=seq.count('Y')))
                reaction.participants.add(val.species_coefficients.get_or_create(coefficient=seq.count('V')))

    def gen_rate_laws(self):

        compartment = self.model.compartments.get_one(id='c')
        submodel = self.submodel

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
                                          k_cat=10,
                                          k_m=1)
