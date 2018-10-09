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

from scipy.constants import Avogadro
import numpy as np
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
            rxn.participants.add(monopp.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient = -100))
            rxn.participants.add(monopp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient =  100))

            #Create conversion reactions
            rxn = submodel.reactions.get_or_create(id='conversion_'+monopp.id+'_'+tripp.id)
            rxn.participants = []
            rxn.participants.add(monopp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient=-100))
            rxn.participants.add(tripp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient = 100))

        # Generate reactions associated with tRna/AA maintenece
        for observable_kb in cell.observables:
            if observable_kb.id[0:5] =='tRNA_':
                rxn = submodel.reactions.get_or_create(id='transfer_'+observable_kb.name)
                rxn.participants = []

                for tRNA_specie_kb in observable_kb.species:
                    tRNA_specie_type_model = model.species_types.get_one(id=tRNA_specie_kb.species.species_type.id)
                    tRNA_model = tRNA_specie_type_model.species.get_one(compartment=c)

                    rxn.participants.add(
                        tRNA_model.species_type.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient = -100))
                    rxn.participants.add(
                        tRNA_model.species_type.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient = 100))

        # Generate reactions associated with H maintenece
        rxn = submodel.reactions.get_or_create(id='transfer_h')
        rxn.participants = []
        rxn.participants.add(h_type.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient=-100))
        rxn.participants.add(h_type.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient= 100))

    def gen_rate_laws(self):
        """ Generate rate laws associated with min metabolism model

            Raises:
                :obj:`ValueError:` if any phosphate species are missing from the model
        """
        # If need X reactions / s, then rate = X/(V*N_Avogadro)
        # TODO: change flat rate to match xTP/AA consumption

        cell = self.knowledge_base.cell
        model = self.model
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')
        volume = cell.properties.get_one(id='initial_volume').value
        cc_length = cell.properties.get_one(id='cell_cycle_length').value

        monopp_trasfer_rate = self.calc_monopp_trasfer_rate()
        H_trasfer_rate = self.calc_H_trasfer_rate()
        monopp_conversion_rate = self.calc_monopp_conversion_rate()

        for rxn in submodel.reactions:

            # Reactions imported from KB have rate laws defined
            if rxn.id[-7:]=='_fromKB':
                continue

            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward

            # Rates for transfer reactions form extracellular space
            if rxn.id[0:9]=='transfer_':

                if rxn.id[0:13]=='transfer_tRNA':
                    expression='0.0000000000589232'
                    if 'transfer_tRNA_rate_equation' not in locals():
                        transfer_tRNA_rate_equation = wc_lang.RateLawEquation(expression=expression)
                    rate_law.equation = transfer_tRNA_rate_equation

                elif rxn.id[-2:]=='mp':
                    expression = str(monopp_trasfer_rate) #'0.0000000007852'
                    if 'transfer_xMP_rate_equation' not in locals():
                        transfer_xMP_rate_equation = wc_lang.RateLawEquation(expression=expression)
                    rate_law.equation = transfer_xMP_rate_equation

                elif rxn.id[0:10]=='transfer_h':
                    expression= str(H_trasfer_rate) #'0.000000013058175578434628529'
                    transfer_H_rate_equation = wc_lang.RateLawEquation(expression=expression)
                    rate_law.equation = transfer_H_rate_equation

                else:
                    raise Exception('Invalid transfer reaction id, no associated rate law.')

            elif rxn.id[0:11]=='conversion_':
                expression= str(monopp_conversion_rate) #'0.000000001264589'
                if 'conversion_rate_law_equation' not in locals():
                    conversion_rate_law_equation = wc_lang.RateLawEquation(expression=expression)
                rate_law.equation = conversion_rate_law_equation

            else:
                raise Exception('Invalid reaction id, no associated rate law.')

    def calc_monopp_trasfer_rate(self):
        """ Calculates the rate of monophosphate transfer from the extracellular space """

        cell = self.knowledge_base.cell
        model = self.model

        cc_length = cell.properties.get_one(id='cell_cycle_length').value
        volume = cell.properties.get_one(id='initial_volume').value
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')

        # First calculate the # of TP molecules tied up in RNAs, which is equal to the # of new TP molecules needed
        # to duplicate over the cell cycle.
        # TODO: Currently volume is treated as equal through CC

        n_TP = 0
        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rna_kbs:
            n_TP += rna_kb.get_seq().count('A')
            n_TP += rna_kb.get_seq().count('C')
            n_TP += rna_kb.get_seq().count('G')
            n_TP += rna_kb.get_seq().count('U')

        # Number of monopp transfer reactions needed per sec
        n_rxn_per_sec = (n_TP/4)/cc_length

        # If need X reactions / s, then rate = X/(V*N_Avogadro)
        # Each transfer reactions transports 100 TPs, thus the final /100
        mpp_transfer_rate = (n_rxn_per_sec/(Avogadro * volume))/100

        print('Monopp transfer rate is: ', mpp_transfer_rate, '(7.85e-10)')
        return mpp_transfer_rate

    def calc_H_trasfer_rate(self):
        """ Calculates the rate of H transfer from the extracellular space """

        cell = self.knowledge_base.cell
        model = self.model

        rnas = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        cc_length = cell.properties.get_one(id='cell_cycle_length').value
        volume = cell.properties.get_one(id='initial_volume').value
        submodel = model.submodels.get_one(id='transcription')
        cytosol_kb = cell.compartments.get_one(id='c')

        # First calculate how many H are needed on avg per transcription reaction
        h_per_transcription=[]
        for rxn in submodel.reactions:
            for part in rxn.participants:
                if part.species.species_type.id=='h':
                    h_per_transcription.append(abs(part.coefficient))

        avg_h_per_transcription = np.mean(h_per_transcription)

        # Calculate the average number of transcription reactions needed for cell, same as the
        # number of initial RNAs
        rna_copy_numbers=[]
        for rna in rnas:
            conc = rna.species.get_one(compartment = cytosol_kb).concentrations.value
            rna_copy_numbers.append(round(conc*volume*Avogadro))

        n_transcription_rxns = np.sum(rna_copy_numbers)

        # Calculate the number of H molecules needed to be transfered over the CC
        n_h_transfer = n_transcription_rxns*avg_h_per_transcription
        print('n_h_transfer: ', n_h_transfer)

        # Number of transfer reactions needed per sec
        n_rxn_per_sec = n_h_transfer/cc_length

        # If need X reactions / s, then rate = X/(V*N_Avogadro)
        # Each transfer reactions transports 100 TPs, thus the final /100
        h_transfer_rate = (n_rxn_per_sec/(volume*Avogadro))/100

        print('H transfer rate is: ', h_transfer_rate, '(1.31e-8)')
        return h_transfer_rate

    def calc_monopp_conversion_rate(self):
        """ Calculates the rate of conversion from mono- to triphosphate molecules """

        cell = self.knowledge_base.cell
        model = self.model

        rnas = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        cc_length = cell.properties.get_one(id='cell_cycle_length').value
        volume = cell.properties.get_one(id='initial_volume').value
        submodel = model.submodels.get_one(id='rna_degradation')
        cytosol_kb = cell.compartments.get_one(id='c')

        # Calculate the number of degradation reactions over CC
        n_deg_rxns=0
        for rna in rnas:
            conc = rna.species.get_one(compartment = cytosol_kb).concentrations.value
            rna_copy_number = round(conc*volume*Avogadro)
            n_deg_rxns += ((cc_length/rna.half_life)*rna_copy_number)

        # Calculate the average number of monopps produced per degradation reaction
        mpp_per_deg_rxn = []
        for rxn in submodel.reactions:
            for part in rxn.participants:
                if part.species.species_type.id=='amp':
                    mpp_per_deg_rxn.append(part.coefficient)
                elif part.species.species_type.id=='cmp':
                    mpp_per_deg_rxn.append(part.coefficient)
                elif part.species.species_type.id=='gmp':
                    mpp_per_deg_rxn.append(part.coefficient)
                elif part.species.species_type.id=='ump':
                    mpp_per_deg_rxn.append(part.coefficient)

        avg_mpp_per_deg_rxn = (np.mean(mpp_per_deg_rxn)/4)/len(submodel.reactions)

        # Calculate the number of mpps needed to be converted
        n_mpp_to_convert = avg_mpp_per_deg_rxn*n_deg_rxns

        # Calculate the number of mpps to convert per sec
        n_mpp_to_convert_per_sec = n_mpp_to_convert/cc_length

        # If need X reactions / s, then rate = X/(V*N_Avogadro)
        # Each transfer reactions transports 100 TPs, thus the final /100
        mpp_to_convert_rate = (n_mpp_to_convert_per_sec/(volume*Avogadro))/100

        # Finally, the mpps from the external medium should be converted as well
        mpp_conversion_rate = mpp_to_convert_rate + self.calc_monopp_trasfer_rate()

        print('PP conversion rate is: ', mpp_conversion_rate, '(1.26E-09)')
        return mpp_conversion_rate
