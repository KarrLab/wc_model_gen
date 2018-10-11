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

        mpps = {}
        for id in ['amp', 'cmp', 'gmp', 'ump']:
            mpps[id] = model.species_types.get_one(id=id)

        # Confirm that all phosphate species were found
        missing = []
        for d in [tripps, mpps]:
            for id, pps_type in d.items():
                if pps_type is None:
                    missing.append(id)
        if missing:
            raise ValueError("'{}' not found in model.species".format(', '.join(missing)))

        # Generate reactions associated with nucleophosphate maintenance
        for mpp, tripp in zip(mpps.values(), tripps.values()):

            # Create transfer reaction
            rxn = submodel.reactions.get_or_create(id='transfer_'+mpp.id)
            rxn.participants = []
            rxn.participants.add(mpp.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient = -10))
            rxn.participants.add(mpp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient =  10))

            #Create conversion reactions
            rxn = submodel.reactions.get_or_create(id='conversion_'+mpp.id+'_'+tripp.id)
            rxn.participants = []
            rxn.participants.add(mpp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient=-10))
            rxn.participants.add(tripp.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient = 10))

        # Generate reactions associated with tRna/AA maintenece
        for observable_kb in cell.observables:
            if observable_kb.id[0:5] =='tRNA_':
                rxn = submodel.reactions.get_or_create(id='transfer_'+observable_kb.name)
                rxn.participants = []

                for tRNA_specie_kb in observable_kb.species:
                    tRNA_specie_type_model = model.species_types.get_one(id=tRNA_specie_kb.species.species_type.id)
                    tRNA_model = tRNA_specie_type_model.species.get_one(compartment=c)

                    rxn.participants.add(
                        tRNA_model.species_type.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient = -10))
                    rxn.participants.add(
                        tRNA_model.species_type.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient = 10))

        # Generate reactions associated with H maintenece
        rxn = submodel.reactions.get_or_create(id='transfer_h')
        rxn.participants = []
        rxn.participants.add(h_type.species.get_or_create(compartment=e).species_coefficients.get_or_create(coefficient=-10))
        rxn.participants.add(h_type.species.get_or_create(compartment=c).species_coefficients.get_or_create(coefficient= 10))

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

        mpp_trasfer_rate = self.calc_mpp_trasfer_rate()
        H_trasfer_rate = self.calc_H_trasfer_rate()
        mpp_conversion_rate = self.calc_mpp_conversion_rate() + mpp_trasfer_rate

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
                    expression = str(mpp_trasfer_rate) #'0.0000000007852'
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
                expression= str(mpp_conversion_rate) #'0.000000001264589'
                if 'conversion_rate_law_equation' not in locals():
                    conversion_rate_law_equation = wc_lang.RateLawEquation(expression=expression)
                rate_law.equation = conversion_rate_law_equation

            else:
                raise Exception('Invalid reaction id, no associated rate law.')

    def calc_H_trasfer_rate(self):
        """ Calculates the rate of H transfer from the extracellular space """

        # Calculate # of H needed on avg per transcription reaction
        avg_H_per_transcription = self.calc_avg_H_per_transcription()

        # Calculate avg # of transcription reactions over CC = # of initial RNAs + # of degrad rxns
        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        n_degrad_rxns = self.calc_expected_n_degrad_rxns()
        n_rnas = self.calc_avg_rna_copy_num()*len(rnas_kb)
        n_transcription_rxns = n_rnas + n_degrad_rxns

        # Calculate the # of H molecules needed to be transfered over the CC
        n_h_transfer = (n_transcription_rxns-n_degrad_rxns)*avg_H_per_transcription

        # Each transfer reactions transports 10 TPs, thus the final /10
        cc_length = self.knowledge_base.cell.properties.get_one(id='cell_cycle_length').value
        h_transfer_rate = n_h_transfer/cc_length/10
        return h_transfer_rate

    def calc_mpp_trasfer_rate(self):
        """ Calculates the rate of monophosphate transfer from the extracellular space """

        # Calculate # of new TPP molecules needed = # of TPP molecules tied up in RNAs at T=0
        avg_tpp_per_rna = self.calc_avg_tpp_per_rna()
        avg_rna_copy_num = self.calc_avg_rna_copy_num()

        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        tpp_in_cell = avg_tpp_per_rna*avg_rna_copy_num*len(rnas_kb)

        # Each transfer reactions transports 10 TPs, thus the final /10
        cc_length = self.knowledge_base.cell.properties.get_one(id='cell_cycle_length').value
        mpp_transfer_rate = (tpp_in_cell/cc_length/10)*1.075 # Slight boost added due to low initial conc
        return mpp_transfer_rate

    def calc_mpp_conversion_rate(self):
        """ Calculates the rate of conversion from mono- to triphosphate molecules """

        # Calculate expected # of degradation reactions over CC
        n_deg_rxns  = self.calc_expected_n_degrad_rxns()

        # Calculate avg # of mpps produced per degrad rxns; same as TPP content of RNAs
        avg_mpp_per_deg_rxn = self.calc_avg_tpp_per_rna()

        # Calculate # of mpps needed to be converted
        n_mpp_to_convert = avg_mpp_per_deg_rxn*n_deg_rxns

        # Each transfer reactions transports 10 TPs, thus the final /10
        cc_length = self.knowledge_base.cell.properties.get_one(id='cell_cycle_length').value
        mpp_conversion_rate = n_mpp_to_convert/cc_length/10
        return mpp_conversion_rate # This is only the rate from degradation, needs to add new mpp conversion!

    """ auxiliary functions """
    def calc_avg_tpp_per_rna(self):
        """ Calculates the average triphosphate content of RNAs """

        cell = self.knowledge_base.cell
        model = self.model

        n_tpp = 0
        rna_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rna_kbs:
            n_tpp += rna_kb.get_seq().count('A')
            n_tpp += rna_kb.get_seq().count('C')
            n_tpp += rna_kb.get_seq().count('G')
            n_tpp += rna_kb.get_seq().count('U')

        avg_tpp_per_rna = n_tpp/4/len(rna_kbs)
        return avg_tpp_per_rna

    def calc_avg_H_per_transcription(self):
        """ Calculates the average H needed for a transcription reaction """
        submodel = self.model.submodels.get_one(id='transcription')

        h_per_transcription=[]
        for rxn in submodel.reactions:
            for part in rxn.participants:
                if part.species.species_type.id=='h':
                    h_per_transcription.append(abs(part.coefficient))

        avg_H_per_transcription = np.mean(h_per_transcription)
        return avg_H_per_transcription

    def calc_expected_n_degrad_rxns(self):
        """ Calculates the expected # of RNA degradation reactions over the CC"""

        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        cc_length = self.knowledge_base.cell.properties.get_one(id='cell_cycle_length').value
        volume = self.knowledge_base.cell.properties.get_one(id='initial_volume').value
        cytosol_kb = self.knowledge_base.cell.compartments.get_one(id='c')

        # Calculate the # of degradation reactions over CC
        # Each Rna degrad reaction fires: (cc_length/rna.half_life)*rna_copy_number(t=0)
        n_deg_rxns=0
        for rna in rnas_kb:
            conc = rna.species.get_one(compartment = cytosol_kb).concentrations.value
            rna_copy_num = round(conc*volume*Avogadro)
            n_deg_rxns += ((cc_length/rna.half_life)*rna_copy_num)

        return n_deg_rxns

    def calc_avg_rna_copy_num(self):
        """ Calculates the # of RNA molecules at t=0 """

        volume = self.knowledge_base.cell.properties.get_one(id='initial_volume').value
        cytosol_kb = self.knowledge_base.cell.compartments.get_one(id='c')
        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)

        rna_copy_num=[]
        for rna in rnas_kb:
            conc = rna.species.get_one(compartment = cytosol_kb).concentrations.value
            rna_copy_num.append(round(conc*volume*Avogadro))

        avg_rna_copy_num = np.mean(rna_copy_num)
        return avg_rna_copy_num
