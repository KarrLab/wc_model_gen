""" Generator for metabolism submodels based on KBs for random in silico organisms
:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
         Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
         Arthur Goldberg <Arthur.Goldberg@mssm.edu>

:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT

TODO: improve terminology to better distinguish this and the metabolism species generation
"""

import wc_model_gen
import wc_lang
import wc_kb


class MetabolismSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for metabolism submodel """

    def gen_reactions(self):
        """ Generate reactions assocated with min model

        Raises:
            :obj:`ValueError:` if any phosphate species are missing from the model
        """
        cell = self.knowledge_base.cell
        model = self.model
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')
        h_type = model.species_types.get_one(id='h')

        # Get species involved in reaction
        tripps = {}
        for id in ['atp', 'ctp', 'gtp', 'utp']:
            tripps[id] = model.species_types.get_one(id=id)

        monopps = {}
        for id in ['amp', 'cmp', 'gmp', 'ump']:
            monopps[id] = model.species_types.get_one(id=id)

        # Confirm that all phosphate species were found
        missing = []
        for d in [tripps, monopps]:
            for id, pps_type in d.items():
                if pps_type is None:
                    missing.append(id)
        if missing:
            raise ValueError("'{}' not found in model.species".format(', '.join(missing)))

        # Generate reactions associated with nucleophosphate maintenance
        for monopp, tripp in zip(monopps.values(), tripps.values()):

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

                    rxn.participants.add(tRNA_model.species_type.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient=-300))
                    rxn.participants.add(tRNA_model.species_type.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient=300))

        # Generate reactions associated with H maintenece
        rxn = submodel.reactions.get_or_create(id='transfer_h')
        rxn.participants = []
        rxn.participants.add(h_type.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient=-300))
        rxn.participants.add(h_type.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient=300))

    def gen_rate_laws(self):
        """ Generate rate laws associated with min metabolism model """
        # If need X reactions / s, then rate = X/(V*N_Avogadro)
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

                # rate for tRNA transfer reactions
                if rxn.id[0:13]=='transfer_tRNA':
                    expression='0.00000000005'

                    if 'transfer_tRNA_rate_equation' not in locals():
                        transfer_tRNA_rate_equation = wc_lang.RateLawEquation(expression=expression)
                    rate_law.equation = transfer_tRNA_rate_equation

                # rate for H transfer reactions; CALIB
                elif rxn.id[0:10]=='transfer_h':
                    expression='0.000000013058175578434628529'
                    if 'transfer_H_rate_equation' not in locals():
                        transfer_H_rate_equation = wc_lang.RateLawEquation(expression=expression)
                    rate_law.equation = transfer_H_rate_equation

                # rate for xMP molecules
                else:
                    expression='0.0000000007852'
                    if 'transfer_H_rate_law_equation' not in locals():
                        transfer_H_rate_law_equation = wc_lang.RateLawEquation(expression=expression)
                    rate_law.equation = transfer_H_rate_law_equation

            elif rxn.id[0:11]=='conversion_':
                expression='0.000000001'
                #expression='0.0000026'
                if 'conversion_rate_law_equation' not in locals():
                    conversion_rate_law_equation = wc_lang.RateLawEquation(expression=expression)
                rate_law.equation = conversion_rate_law_equation

            else:
                raise Exception('Unknown reaction type')
