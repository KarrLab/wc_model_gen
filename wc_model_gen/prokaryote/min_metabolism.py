""" Generator for metabolism submodels based on KBs for random in silico organisms
:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
         Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>

:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT

TODO: improve temrinology to better distinguish this and the metabolism species generation
"""

import wc_model_gen
import wc_lang
import wc_kb


class MinMetabolismSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for metabolism submodel """

    def gen_compartments(self):
        pass

    def gen_parameters(self):
        pass

    def gen_species(self):
        pass

    def gen_reactions(self):
        """ Generate reactions assocated with min model """
        cell = self.knowledge_base.cell
        model = self.model
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')

        # Get species involved in reaction
        atp_type = model.species_types.get_one(id='atp')
        ctp_type = model.species_types.get_one(id='ctp')
        gtp_type = model.species_types.get_one(id='gtp')
        utp_type = model.species_types.get_one(id='utp')
        amp_type = model.species_types.get_one(id='amp')
        cmp_type = model.species_types.get_one(id='cmp')
        gmp_type = model.species_types.get_one(id='gmp')
        ump_type = model.species_types.get_one(id='ump')

        #.species.get_one(compartment=cytosol)
        tripps  = [atp_type, ctp_type, gtp_type ,utp_type]
        monopps = [amp_type, cmp_type, gmp_type, ump_type]

        # Generate reactions associated with nucleophosphate maintenece
        for monopp, tripp in zip(monopps, tripps):

            # Create transfer reaction
            rxn = submodel.reactions.get_or_create(id='transfer_'+monopp.id)
            rxn.participants = []
            rxn.participants.add(monopp.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(monopp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient=1))

            #Create conversion reactions
            rxn = submodel.reactions.get_or_create(id='conversion_'+monopp.id+'_'+tripp.id)
            rxn.participants = []
            rxn.participants.add(monopp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(tripp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient=1))

        # Generate reactions associated with tRna/AA maintenece
        for observable_kb in cell.observables:
            if observable_kb.id[0:5] =='tRNA_':
                rxn = submodel.reactions.get_or_create(id='transfer_'+observable_kb.name)
                rxn.participants = []

                for tRNA_specie_kb in observable_kb.species:
                    tRNA_specie_type_model = model.species_types.get_one(id=tRNA_specie_kb.species.species_type.id)
                    tRNA_model = tRNA_specie_type_model.species.get_one(compartment=c)

                    rxn.participants.add(tRNA_model.species_type.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient=-1))
                    rxn.participants.add(tRNA_model.species_type.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient=1))

    def gen_rate_laws(self):
        """ Generate rate laws associated with min metabolism model """
        # rate of 0.0000000000004151 ~ 1 reaction /s
        # TODO: change flat rate to match xTP/AA consumption

        cell = self.knowledge_base.cell
        model = self.model
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')            

        for rxn in submodel.reactions:
            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward

            if rxn.id[0:9]=='transfer_':
                expression='0.0000000000005'
                if 'rate_law_equation' not in locals():
                    rate_law_equation = wc_lang.RateLawEquation(expression=expression)

                rate_law.equation = rate_law_equation

            elif rxn.id[0:11]=='conversion_':
                expression='0.0000000000006'
                if 'rate_law_equation' not in locals():
                    rate_law_equation = wc_lang.RateLawEquation(expression=expression)

                rate_law.equation = rate_law_equation

            else:
                raise Exception('Unknown reaction type')
