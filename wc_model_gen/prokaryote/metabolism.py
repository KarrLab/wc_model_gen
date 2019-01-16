""" Generator for metabolism submodels based on KBs for random in silico organisms

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <karr@mssm.edu>
:Author: Ashwin Srinivasan <ashwins@mit.edu>
:Author: Arthur Goldberg <Arthur.Goldberg@mssm.edu>
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT

TODO: improve terminology to better distinguish this and the metabolism species generation
"""

from scipy.constants import Avogadro
from wc_utils.util.ontology import wcm_ontology
import numpy
import wc_lang
import wc_kb
import wc_model_gen


class MetabolismSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for metabolism submodel """

    # self.reaction_scale is the number of molecules transfered / converted with each reaction.
    # Increasing the number reduces the number of metabolic reactions, significantly reducing simulation time
    reaction_scale = 100

    def gen_reactions(self):
        """ Generate reactions assocated with min model

        Raises:
            :obj:`ValueError:` if any phosphate species are missing from the model
        """

        # Get species, model components
        cell = self.knowledge_base.cell
        model = self.model
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')
        h_type = model.species_types.get_one(id='h')
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=c)
        h = model.species_types.get_one(id='h').species.get_one(compartment=c)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=c)
        adp = model.species_types.get_one(id='adp').species.get_one(compartment=c)

        # Get species involved in reaction
        mpps = {}
        for id in ['ump', 'amp', 'cmp', 'gmp']:
            mpps[id] = model.species_types.get_one(id=id)

        dpps = {}
        for id in ['udp', 'adp', 'cdp', 'gdp']:
            dpps[id] = model.species_types.get_one(id=id)

        tpps = {}
        for id in ['utp', 'atp', 'ctp', 'gtp']:
            tpps[id] = model.species_types.get_one(id=id)

        # Confirm that all phosphate species were found
        missing = []
        for d in [tpps, dpps, mpps]:
            for id, pps_type in d.items():
                if pps_type is None:
                    missing.append(id)
        if missing:
            raise ValueError("'{}' not found in model.species".format(', '.join(missing)))

        # Generate reactions associated with nucleophosphate maintenance
        for dpp, tripp in zip(dpps.values(), tpps.values()):

            # get/create species
            dpp_e   = model.species.get_or_create(species_type=dpp,   compartment=e)
            dpp_c   = model.species.get_or_create(species_type=dpp,   compartment=c)
            tripp_c = model.species.get_or_create(species_type=tripp, compartment=c)
            dpp_e.id = dpp_e.gen_id()
            dpp_c.id = dpp_c.gen_id()
            tripp_c.id = tripp_c.gen_id()

            # Create transfer reaction
            rxn = model.reactions.get_or_create(submodel=submodel, id='transfer_'+dpp.id)
            rxn.participants = []
            rxn.participants.add(dpp_e.species_coefficients.get_or_create(coefficient=-self.reaction_scale))
            rxn.participants.add(dpp_c.species_coefficients.get_or_create(coefficient=self.reaction_scale))

            # Create conversion reactions
            rxn = model.reactions.get_or_create(submodel=submodel, id='conversion_'+dpp.id+'_'+tripp.id)
            rxn.participants = []
            rxn.participants.add(dpp_c.species_coefficients.get_or_create(coefficient=-self.reaction_scale))
            rxn.participants.add(pi.species_coefficients.get_or_create(coefficient=-self.reaction_scale))
            #rxn.participants.add(h.species_coefficients.get_or_create(coefficient=-self.reaction_scale))

            rxn.participants.add(tripp_c.species_coefficients.get_or_create(coefficient=self.reaction_scale))
            #rxn.participants.add(h.species_coefficients.get_or_create(coefficient=self.reaction_scale))
            rxn.participants.add(h2o.species_coefficients.get_or_create(coefficient=self.reaction_scale))

        # Generate reactions associated with tRna/AA maintenece
        for observable_kb in cell.observables:
            if observable_kb.id[0:5] == 'tRNA_':
                rxn = model.reactions.get_or_create(submodel=submodel, id='transfer_'+observable_kb.name)
                rxn.participants = []

                for tRNA_specie_kb in observable_kb.species:
                    tRNA_specie_type_model = model.species_types.get_one(id=tRNA_specie_kb.species.species_type.id)
                    tRNA_model = tRNA_specie_type_model.species.get_one(compartment=c)

                    trna_e = model.species.get_or_create(
                        species_type=tRNA_model.species_type,
                        compartment=e)
                    trna_c = model.species.get_or_create(
                        species_type=tRNA_model.species_type,
                        compartment=c)
                    trna_e.id = trna_e.gen_id()
                    trna_c.id = trna_c.gen_id()
                    rxn.participants.add(trna_e.species_coefficients.get_or_create(coefficient=-self.reaction_scale))
                    rxn.participants.add(trna_c.species_coefficients.get_or_create(coefficient=self.reaction_scale))

        # Generate reactions associated with H maintenece
        rxn = model.reactions.get_or_create(submodel=submodel, id='transfer_h')
        rxn.participants = []

        h_e = model.species.get_or_create(species_type=h_type, compartment=e)
        h_c = model.species.get_or_create(species_type=h_type, compartment=c)
        h_e.id = h_e.gen_id()
        h_c.id = h_c.gen_id()
        rxn.participants.add(h_e.species_coefficients.get_or_create(coefficient=-self.reaction_scale))
        rxn.participants.add(h_c.species_coefficients.get_or_create(coefficient=self.reaction_scale))

    def gen_rate_laws(self):
        """ Generate rate laws associated with min metabolism model

            Raises:
                :obj:`ValueError:` if any phosphate species are missing from the model
                :obj:`Exception:` if there is a reaction with unexpected ID
        """
        cell = self.knowledge_base.cell
        model = self.model
        submodels = self.model.submodels
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')
        volume = cell.properties.get_one(id='mean_volume').value
        cc_length = cell.properties.get_one(id='mean_doubling_time').value

        # Calculate rates
        mpp_transfer_rate = self.calc_mpp_transfer_rate()
        mpp_conversion_rate = self.calc_mpp_conversion_rate() + mpp_transfer_rate
        aa_transfer_rate = self.calc_aa_transfer_rate()
        h_transfer_rate = self.calc_h_transfer_rate()

        # Loop through reactions and set rate laws
        for rxn in submodel.reactions:

            # Reactions imported from KB have rate laws defined
            if rxn.id.endswith('_kb'):
                continue

            rate_law = model.rate_laws.create(
                reaction=rxn,
                direction=wc_lang.RateLawDirection.forward)
            rate_law.id = rate_law.gen_id()

            param = model.parameters.create(id='k_cat_{}'.format(rxn.id),
                                            type=wcm_ontology['WCM:k_cat'],
                                            units='s^-1')

            # Rates for transfer reactions from extracellular space
            if rxn.id[0:9] == 'transfer_':
                # tRNA transfer
                if rxn.id[-8:-4] == 'tRNA':
                    param.value = 0.0000000000589232

                # Monophosphate reaction
                elif rxn.id[-2:] == 'dp':
                    param.value = mpp_transfer_rate

                # Hydrogen transfer
                elif rxn.id[-1:] == 'h':
                    param.value = h_transfer_rate

                else:
                    raise Exception('{}: invalid reaction id, no associated rate law is defined.'.format(rxn.id))

            elif rxn.id[0:11] == 'conversion_':
                param.value = mpp_conversion_rate

            else:
                raise Exception('{}: invalid reaction id, no associated rate law is defined.'.format(rxn.id))

            rate_law.expression, error = wc_lang.RateLawExpression.deserialize(param.id, {
                wc_lang.Parameter: {param.id: param},
            })
            assert error is None, str(error)

        # Apply correction terms if translation submodel is present
        #if model.submodels.get_one(id='translation'):
        #    gtp_corr_rate = self.calc_gtp_corr_rate()
        #    model.parameters.get_one(id='k_cat_transfer_gmp').value += gtp_corr_rate
        #    model.parameters.get_one(id='k_cat_conversion_gmp_gtp').value += gtp_corr_rate

    """ Auxiliary functions """

    def calc_h_transfer_rate(self):
        """ Calculates the rate of H transfer from the extracellular space """
        kb = self.knowledge_base

        # Calculate # of H needed on avg per transcription reaction
        avg_H_per_transcription = self.calc_h_per_transcript()

        # Calculate avg # of transcription reactions over CC = # of initial RNAs + # of degrad rxns
        rnas_kb = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        n_degrad_rxns = self.calc_rna_degrad_rxns()
        n_rnas = self.calc_rna_copy_num() * len(rnas_kb)
        n_transcription_rxns = n_rnas + n_degrad_rxns

        # Calculate the # of H molecules needed to be transfered over the CC
        n_h_transfer = (n_transcription_rxns-n_degrad_rxns)*avg_H_per_transcription

        # Each transfer reactions transports 10 TPs, thus the final /10
        cc_length = kb.cell.properties.get_one(id='mean_doubling_time').value
        h_transfer_rate = n_h_transfer / cc_length / self.reaction_scale
        #print('h_transfer_rate: ', h_transfer_rate)
        return h_transfer_rate

    def calc_mpp_transfer_rate(self):
        """ Calculates the rate of monophosphate transfer from the extracellular space """

        # Calculate # of new TPP molecules needed = # of TPP molecules tied up in RNAs at T=0
        avg_tpp_per_rna = self.calc_tpp_per_rna()
        avg_rna_copy_num = self.calc_rna_copy_num()

        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        tpp_in_cell = avg_tpp_per_rna*avg_rna_copy_num*len(rnas_kb)

        cc_length = self.knowledge_base.cell.properties.get_one(id='mean_doubling_time').value
        mpp_transfer_rate = (tpp_in_cell/cc_length/self.reaction_scale)

        #print('mpp_transfer_rate: ', mpp_transfer_rate)
        return mpp_transfer_rate

    def calc_aa_transfer_rate(self):
        """ Calculates the rate of monophosphate transfer from the extracellular space """

        # Calculate # of new AA molecules needed = # of AA molecules tied up in prots at t=0
        avg_aa_per_prot = self.calc_aa_per_prot()
        avg_prot_copy_num = self.calc_prot_copy_num()

        prots_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        aa_in_cell = avg_aa_per_prot*avg_prot_copy_num*len(prots_kb)

        cc_length = self.knowledge_base.cell.properties.get_one(id='mean_doubling_time').value
        aa_transfer_rate = (aa_in_cell/cc_length/self.reaction_scale)

        #print('aa_transfer_rate: ', aa_transfer_rate)
        return aa_transfer_rate

    def calc_mpp_conversion_rate(self):
        """ Calculates the rate of conversion from mono- to triphosphate molecules """
        # Calculate expected # of degradation reactions over CC
        n_rna_deg_rxns = self.calc_rna_degrad_rxns()

        # Calculate avg # of mpps produced per degrad rxns; same as TPP content of RNAs
        avg_mpp_per_deg_rxn = self.calc_tpp_per_rna()

        # Calculate # of mpps needed to be converted
        n_mpp_to_convert = avg_mpp_per_deg_rxn*n_rna_deg_rxns

        cc_length = self.knowledge_base.cell.properties.get_one(id='mean_doubling_time').value
        mpp_conversion_rate = n_mpp_to_convert/cc_length/self.reaction_scale

        #print('mpp_conversion_rate: ', mpp_conversion_rate)
        return mpp_conversion_rate  # This is only the rate from degradation, needs to add new mpp conversion!

    def calc_tpp_per_rna(self):
        """ Calculates the average triphosphate content of RNAs """

        cell = self.knowledge_base.cell
        model = self.model

        n_tpp = 0
        rnas_kb = cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)
        for rna_kb in rnas_kb:
            n_tpp += rna_kb.get_seq().count('A')
            n_tpp += rna_kb.get_seq().count('C')
            n_tpp += rna_kb.get_seq().count('G')
            n_tpp += rna_kb.get_seq().count('U')

        avg_tpp_per_rna = n_tpp/4/len(rnas_kb)

        #print('avg_tpp_per_rna: ', avg_tpp_per_rna)
        return avg_tpp_per_rna

    def calc_aa_per_prot(self):
        """ Calculates the average amino acid content of proteins """
        cell = self.knowledge_base.cell
        model = self.model

        proteins_kb = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        n_aa = 0
        for protein_kb in proteins_kb:
            n_aa += len(protein_kb.get_seq())

        observables_kb = cell.observables
        n_trnas = 0  # we need to count n of TRNAs, since in reduced models there are only 4
        for observable in observables_kb:
            if observable.id[0:5] == 'tRNA_':
                n_trnas += 1

        avg_aa_per_prot = n_aa/(n_trnas*len(proteins_kb))

        #print('avg_aa_per_prot: ', avg_aa_per_prot)
        return avg_aa_per_prot

    def calc_h_per_transcript(self):
        """ Calculates the average H needed for a transcription reaction """
        submodel = self.model.submodels.get_one(id='transcription')

        h_per_transcription = []
        for rxn in submodel.reactions:
            for part in rxn.participants:
                if part.species.species_type.id == 'h':
                    h_per_transcription.append(abs(part.coefficient))

        avg_H_per_transcription = numpy.mean(h_per_transcription)

        #print('avg_H_per_transcription: ', avg_H_per_transcription)
        return avg_H_per_transcription

    def calc_gtp_per_translate(self):
        """ Calculates the average GTP needed for a translation reaction """
        submodel = self.model.submodels.get_one(id='translation')

        gtp_per_translation = []
        for rxn in submodel.reactions:
            for part in rxn.participants:
                if part.species.species_type.id == 'gtp':
                    gtp_per_translation.append(abs(part.coefficient))

        avg_gtp_per_translate = numpy.mean(gtp_per_translation)

        #print('avg_gtp_per_translate: ', avg_gtp_per_translate)
        return avg_gtp_per_translate

    def calc_rna_degrad_rxns(self):
        """ Calculates the expected # of RNA degradation reactions over the CC"""

        cytosol_lang = self.model.compartments.get_one(id='c')
        cytosol_kb = self.knowledge_base.cell.compartments.get_one(id='c')
        rnas_lang = self.model.species_types.get(type=wcm_ontology['WCM:RNA']) # RNA
        cc_length = self.knowledge_base.cell.properties.get_one(id='mean_doubling_time').value
        volume = self.knowledge_base.cell.properties.get_one(id='mean_volume').value

        # Calculate the # of degradation reactions over CC
        # Each Rna degrad reaction fires: (cc_length/rna.half_life)*rna_copy_number(t=0)
        n_rna_deg_rxns = 0
        for rna_lang in rnas_lang:
            rna_kb = self.knowledge_base.cell.species_types.get_one(id=rna_lang.id)
            half_life = rna_kb.half_life

            conc = rna_lang.species.get_one(compartment=cytosol_lang) \
                .distribution_init_concentration.mean
            rna_copy_num = round(conc * volume * Avogadro)
            n_rna_deg_rxns += ((cc_length / half_life) * rna_copy_num)

        #print('n_rna_deg_rxns: ', n_rna_deg_rxns)
        return n_rna_deg_rxns

    def calc_prot_degrad_rxns(self):
        """ Calculates the expected # of RNA degradation reactions over the CC"""

        cytosol_lang = self.model.compartments.get_one(id='c')
        cytosol_kb = self.knowledge_base.cell.compartments.get_one(id='c')
        prots_lang = self.model.species_types.get(type=wcm_ontology['WCM:protein']) # protein
        cc_length = self.knowledge_base.cell.properties.get_one(id='mean_doubling_time').value
        volume = self.knowledge_base.cell.properties.get_one(id='mean_volume').value

        # Calculate the # of degradation reactions over CC
        # Each prot degrad reaction fires: (cc_length/prot.half_life)*prot_copy_number(t=0)
        n_prot_deg_rxns = 0
        for prot_lang in prots_lang:
            prot_kb = self.knowledge_base.cell.species_types.get_one(id=prot_lang.id)
            if isinstance(prot_kb, wc_kb.core.ComplexSpeciesType):
                continue
            half_life = prot_kb.half_life

            conc = prot_lang.species.get_one(compartment=cytosol_lang)\
                .distribution_init_concentration.mean
            prot_copy_num = round(conc * volume * Avogadro)
            n_prot_deg_rxns += ((cc_length / half_life) * prot_copy_num)

        #n_prot_deg_rxns = round(n_prot_deg_rxns)
        #print('n_prot_deg_rxns: ', n_prot_deg_rxns)
        return n_prot_deg_rxns

    def calc_rna_copy_num(self):
        """ Calculates the # of RNA molecules at t=0 """

        volume = self.knowledge_base.cell.properties.get_one(id='mean_volume').value
        cytosol_kb = self.knowledge_base.cell.compartments.get_one(id='c')
        rnas_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)

        rna_copy_num = []
        for rna in rnas_kb:
            conc = rna.species.get_one(compartment=cytosol_kb).concentration.value
            rna_copy_num.append(round(conc * volume * Avogadro))

        avg_rna_copy_num = numpy.mean(rna_copy_num)

        #print('avg_rna_copy_num: ', avg_rna_copy_num)
        return avg_rna_copy_num

    def calc_prot_copy_num(self):
        """ Calculates the # of RNA molecules at t=0 """
        # TODO: check if conc is in CN already

        volume = self.knowledge_base.cell.properties.get_one(id='mean_volume').value
        cytosol_kb = self.knowledge_base.cell.compartments.get_one(id='c')
        prots_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)

        prot_copy_num = []
        for prot in prots_kb:
            conc = prot.species.get_one(compartment=cytosol_kb).concentration.value
            prot_copy_num.append(round(conc*volume*Avogadro))

        avg_prot_copy_num = numpy.mean(prot_copy_num)

        #print('avg_prot_copy_num: ', avg_prot_copy_num)
        return avg_prot_copy_num

    def calc_gtp_corr_rate(self):
        """ Calculates the extra amount of GTP needed if model has translation subunit """
        cell = self.knowledge_base.cell
        cc_length = cell.properties.get_one(id='mean_doubling_time').value

        avg_gtp_per_translate = self.calc_gtp_per_translate()
        avg_prot_copy_num = self.calc_prot_copy_num()  # Number of new proteins needed
        n_prot_degrad_rxns = self.calc_prot_degrad_rxns()
        total_translation_gtp = avg_gtp_per_translate * (avg_prot_copy_num + n_prot_degrad_rxns)

        gtp_corr_rate = (total_translation_gtp/cc_length/self.reaction_scale)

        #print('gtp_corr_rate: ', gtp_corr_rate)
        return gtp_corr_rate
