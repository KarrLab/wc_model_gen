""" Generator for translation, protein folding and translocation submodel for eukaryotes
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-06-14
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.units import unit_registry
import wc_model_gen.global_vars as gvar
import wc_model_gen.utils as utils
import Bio.Seq
import math
import numpy
import scipy.constants
import wc_kb
import wc_lang
import wc_model_gen
import wc_ontology


class TranslationTranslocationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for translation, protein folding and translocation submodel

        Translation, protein folding and translocation processes are 
        modeled as three reaction steps in this submodel:

        1. Translation initiation where ribosomes and methionine bind to the mRNA. 
           For nuclear mRNAs, transport from the nucleus to the cytoplasm are lumped with
           this reaction. The energetic of met-tRNA charging is included;
        2. Translation elongation and termination are lumped into one reaction that produces
           nascent polypeptides. The energetic of amino-acid-tRNA charging is included;
        3. Protein folding and translocation to each organelle/compartment are lumped into 
           one reaction

        Options:
        * cytoplasmic_ribosome (:obj:`str`): name of cytoplasmic ribosome
        * mitochondrial_ribosome (:obj:`str`): name of mitochondrial ribosome
        * cytoplasmic_initiation_factors (:obj:`list` of :obj:`list`): list of lists of the name of
            initiation factors in the cytoplasm, grouped based on similar functions or classes, 
            the default is an empty list
        * mitochondrial_initiation_factors (:obj:`list` of :obj:`list`): list of lists of the name of
            initiation factors in the mitochondria, grouped based on similar functions or classes, 
            the default is an empty list
        * cytoplasmic_elongation_factors (:obj:`list` of :obj:`list`): list of lists of the name of
            elongation factors in the cytoplasm, grouped based on similar functions or classes, 
            the default is an empty list
        * mitochondrial_elongation_factors (:obj:`list` of :obj:`list`): list of lists of the name of
            elongation factors in the mitochondria, grouped based on similar functions or classes, 
            the default is an empty list
        * cytoplasmic_chaperones (:obj:`list` of :obj:`list`): list of lists of the name of
            chaperones in the cytoplasm, grouped based on similar functions or classes, 
            the default is an empty list
        * mitochondrial_chaperones (:obj:`list` of :obj:`list`): list of lists of the name of
            chaperones in the mitochondria, grouped based on similar functions or classes, 
            the default is an empty list
        * er_chaperones (:obj:`list` of :obj:`list`): list of lists of the name of
            chaperones in the endoplasmic reticulum, grouped based on similar functions or classes, 
            the default is an empty list                
        * amino_acid_id_conversion (:obj:`dict`): a dictionary with amino acid standard ids
            as keys and amino acid metabolite ids as values     
        * codon_table (:obj:`dict`, optional): a dictionary with protein id as key and 
            NCBI identifier for translation table as value, the default is 1 (standard table) 
            for all protein
        * cds (:obj:`bool`, optional): True indicates the sequences of protein are complete CDS,
            the default is True        
        * beta (:obj:`float`, optional): ratio of Michaelis-Menten constant to substrate 
            concentration (Km/[S]) for use when estimating Km values, the default value is 1
        * polysome_fraction (:obj:`dict`, optional): a dictionary with mRNA ids as keys and
            fraction of total cellular ribosomes the mRNA is bound to            
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        if 'cytoplasmic_ribosome' not in options:
            raise ValueError('The name of cytoplasmic ribosome has not been provided')
        else:    
            cytoplasmic_ribosome = options['cytoplasmic_ribosome']

        if 'mitochondrial_ribosome' not in options:
            raise ValueError('The name of mitochondrial ribosome has not been provided')
        else:    
            mitochondrial_ribosome = options['mitochondrial_ribosome']

        cytoplasmic_initiation_factors = options.get('cytoplasmic_initiation_factors', [])
        options['cytoplasmic_initiation_factors'] = cytoplasmic_initiation_factors

        mitochondrial_initiation_factors = options.get('mitochondrial_initiation_factors', [])
        options['mitochondrial_initiation_factors'] = mitochondrial_initiation_factors

        cytoplasmic_elongation_factors = options.get('cytoplasmic_elongation_factors', [])
        options['cytoplasmic_elongation_factors'] = cytoplasmic_elongation_factors

        mitochondrial_elongation_factors = options.get('mitochondrial_elongation_factors', [])
        options['mitochondrial_elongation_factors'] = mitochondrial_elongation_factors

        cytoplasmic_chaperones = options.get('cytoplasmic_chaperones', [])
        options['cytoplasmic_chaperones'] = cytoplasmic_chaperones

        mitochondrial_chaperones = options.get('mitochondrial_chaperones', [])
        options['mitochondrial_chaperones'] = mitochondrial_chaperones

        er_chaperones = options.get('er_chaperones', [])
        options['er_chaperones'] = er_chaperones  

        if 'amino_acid_id_conversion' not in options:
            raise ValueError('The dictionary amino_acid_id_conversion has not been provided')
        else:    
            amino_acid_id_conversion = options['amino_acid_id_conversion']

        codon_table = options.get('codon_table', 1)
        options['codon_table'] = codon_table

        cds = options.get('cds', True)
        options['cds'] = cds     

        beta = options.get('beta', 1.)
        options['beta'] = beta        

        if 'polysome_fraction' not in options:
            raise ValueError('The dictionary polysome_fraction has not been provided')
        else:    
            polysome_fraction = options['polysome_fraction']

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        
        cytoplasmic_ribosome = self.options.get('cytoplasmic_ribosome')
        mitochondrial_ribosome = self.options.get('mitochondrial_ribosome')
        amino_acid_id_conversion = self.options.get('amino_acid_id_conversion')
        codon_table = self.options['codon_table']
        cds = self.options['cds']        
        
        cytosol = model.compartments.get_one(id='c')
        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')
        peroxisome = model.compartments.get_one(id='x')
        er = model.compartments.get_one(id='r')
                           
        # Get metabolite species involved in reaction
        amino_acid_participants = list(amino_acid_id_conversion.values()) 
        other_metabolite_participants = ['atp', 'adp,', 'amp', 'gtp','gdp','pi', 'h2o', 'h']
        metabolites = {}
        for met in amino_acid_participants + other_metabolite_participants:
            met_species_type = model.species_types.get_one(id=met)
            metabolites[met] = {
                'c': met_species_type.species.get_one(compartment=cytosol),
                'm': met_species_type.species.get_one(compartment=mitochondrion)
                }

        print('Start generating translation and translocation submodel...')
        
        # Create initiation and elongation reactions for each mRNA
        init_el_rxn_no = 0
        trans_rxn_no = 0        
        self._allowable_queue_len = {}   
        self._translocation_reactions = {}

        mrna_kbs = [i for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType) \
            if i.type==wc_kb.eukaryote.TranscriptType.mRna]
        
        for mrna_kb in mrna_kbs:
            
            mrna_kb_compartment_id = mrna_kb.species[0].compartment.id
            if mrna_kb_compartment_id == 'n':
                mrna_compartment = nucleus
                translation_compartment = cytosol
                ribosome_complex = model.species_types.get_one(
                    name=cytoplasmic_ribosome).species.get_one(compartment=cytosol)
            else:
                mrna_compartment = translation_compartment = mitochondrion
                ribosome_complex = model.species_types.get_one(
                    name=mitochondrial_ribosome).species.get_one(compartment=mitochondrion)            

            # Create initiation reaction
            methionine = model.species_types.get_one(id=amino_acid_id_conversion['M'])            

            if mrna_kb.id in gvar.transcript_ntp_usage:
                mrna_len = gvar.transcript_ntp_usage[mrna_kb.id]['len']
            else:
                seq = mrna_kb.get_seq()
                mrna_len = len(seq)
                ntp_count = gvar.transcript_ntp_usage[mrna_kb.id] = {
                    'A': seq.upper().count('A'),
                    'C': seq.upper().count('C'),
                    'G': seq.upper().count('G'),
                    'U': seq.upper().count('U'),
                    'len': mrna_len
                    }
            
            ribo_binding_site_species = model.species_types.get_one(
                id='{}_ribosome_binding_site'.format(mrna_kb.id)).species[0]
            self._allowable_queue_len[mrna_kb.id] = (ribo_binding_site_species, 
                ribo_binding_site_species.distribution_init_concentration.mean)

            ribo_bound_species_type = model.species_types.get_or_create(
                id='ribo_bound_{}'.format(mrna_kb.id),
                name='Ribosome bound {}'.format(mrna_kb.name),
                type=onto['WC:pseudo_species'],                
                )
            ribo_bound_species_type.structure = wc_lang.ChemicalStructure(
                empirical_formula = ribosome_complex.species_type.structure.empirical_formula +\
                    methionine.structure.empirical_formula,
                molecular_weight = ribosome_complex.species_type.structure.molecular_weight +\
                    methionine.structure.molecular_weight,
                charge = ribosome_complex.species_type.structure.charge +\
                    methionine.structure.charge,
                )
            ribo_bound_species = model.species.get_or_create(
                species_type=ribo_bound_species_type, compartment=translation_compartment)
            ribo_bound_species.id = ribo_bound_species.gen_id()

            conc_model = model.distribution_init_concentrations.create(
                species=ribo_bound_species,
                units=unit_registry.parse_units('molecule'),
                )
            conc_model.id = conc_model.gen_id()

            init_reaction = model.reactions.create(
                submodel=self.submodel, id='translation_initiation_' + mrna_kb.id,
                name='transcription initiation of ' + mrna_kb.name,
                reversible=False, comments='Set to irreversible to model only the net flux')
            
            # Adding participants to LHS 
            # Include 2 GTP hydrolysis and 1 ATP hydrolysis at the initiation factors
            # Include 1 ATP hydrolysis for the charging of tRNA-met            
            init_reaction.participants.append(
                ribosome_complex.species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(
                ribo_binding_site_species.species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(methionine.species.get_one(
                compartment=translation_compartment).species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(metabolites['h2o'][
                translation_compartment.id].species_coefficients.get_or_create(coefficient=-5))
            init_reaction.participants.append(metabolites['atp'][
                translation_compartment.id].species_coefficients.get_or_create(coefficient=-2))
            init_reaction.participants.append(metabolites['gtp'][
                translation_compartment.id].species_coefficients.get_or_create(coefficient=-2))
            
            # Adding participants to RHS          
            init_reaction.participants.append(
                ribo_bound_species.species_coefficients.get_or_create(
                coefficient=1))
            init_reaction.participants.append(metabolites['h'][
                translation_compartment.id].species_coefficients.get_or_create(coefficient=5))
            init_reaction.participants.append(metabolites['amp'][
                translation_compartment.id].species_coefficients.get_or_create(coefficient=1))
            init_reaction.participants.append(metabolites['adp'][
                translation_compartment.id].species_coefficients.get_or_create(coefficient=1))
            init_reaction.participants.append(metabolites['gdp'][
                translation_compartment.id].species_coefficients.get_or_create(coefficient=2))
            init_reaction.participants.append(metabolites['pi'][
                translation_compartment.id].species_coefficients.get_or_create(coefficient=5))
            
            # Create elongation reaction
            protein_model = model.species_types.get_one(id=mrna_kb.protein.id).species.get_or_create(
                compartment=translation_compartment)
            el_reaction = model.reactions.get_or_create(
                submodel=self.submodel, id='translation_elongation_' + mrna_kb.id,
                name='translation elongation of ' + mrna_kb.name,
                reversible=False, comments='Lumped reaction')

            aa_content = {}

            if mrna_kb.protein.id in gvar.protein_aa_usage:                                        
                for aa, aa_id in amino_acid_id_conversion.items():
                    if gvar.protein_aa_usage[mrna_kb.protein.id][aa]:
                        aa_content[aa_id] = gvar.protein_aa_usage[mrna_kb.protein.id][aa]
            else:
                if codon_table == 1:
                    codon_id = 1
                else:
                    codon_id = codon_table[mrna_kb.protein.id]                                        
                protein_seq = ''.join(i for i in mrna_kb.protein.get_seq(table=codon_id, cds=cds) if i!='*')
                for aa in protein_seq:
                    aa_id = amino_acid_id_conversion[aa]
                    if aa_id not in aa_content:
                        aa_content[aa_id] = 1
                    else:
                        aa_content[aa_id] += 1
            
            aa_content[amino_acid_id_conversion['M']] -= 1

            # Adding participants to LHS
            # Include 1 ATP hydrolysis for each tRNA-aa charging
            # Include 1 GTP hydrolysis for each peptide bond formation
            # Include 1 GTP hydrolysis at termination 
            el_reaction.participants.append(ribo_bound_species.species_coefficients.get_or_create(
                coefficient=-1))
            el_reaction.participants.append(metabolites['gtp'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=-len(protein_seq)))
            el_reaction.participants.append(metabolites['atp'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=-(len(protein_seq)-1)))
            el_reaction.participants.append(metabolites['h2o'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=-((len(protein_seq)-1)*2 + 1)))
            for aa_met, count in aa_content.items():                
                el_reaction.participants.append(metabolites[aa_met][
                    translation_compartment.id].species_coefficients.get_or_create(
                        coefficient=-protein_seq.count(aa_std)))            
            
            # Adding participants to RHS
            el_reaction.participants.append(ribosome_complex.species_coefficients.get_or_create(
                coefficient=1))
            el_reaction.participants.append(ribo_binding_site_species.species_coefficients.get_or_create(
                coefficient=1))
            el_reaction.participants.append(protein_model.species_coefficients.get_or_create(
                coefficient=1))
            el_reaction.participants.append(metabolites['amp'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=len(protein_seq)-1))
            el_reaction.participants.append(metabolites['gdp'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=len(protein_seq)))
            el_reaction.participants.append(metabolites['pi'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=(len(protein_seq)-1)*3 + 1))
            el_reaction.participants.append(metabolites['h'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=len(protein_seq)))

            init_el_rxn_no += 1
            
            # Create translocation reactions
            all_localized_comp = [i.compartment for i in model.species_types.get_one(
                id=mrna_kb.protein.id).species if i.compartment!=translation_compartment]
            self._translocation_reactions[mrna_kb] = {}
            for compartment in all_localized_comp:
                
                trans_reaction = model.reactions.get_or_create(
                    submodel=self.submodel, id='translocation_{}_{}_to_{}'.format(
                        mrna_kb.protein.id, translation_compartment.id ,compartment.id),
                    name='translocation of {} from {} to {}'.format(
                        mrna_kb.protein.name, translation_compartment.name, compartment.name),
                    reversible=False, comments='Lumped reaction')
                self._translocation_reactions[mrna_kb][compartment] = trans_reaction

                if compartment.id=='n'
                    energy_compartment = nucleus
                    energy_reactant = 'gtp'
                    energy_product = 'gdp'
                elif compartment.id=='m':
                    energy_compartment = mitochondrion
                    energy_reactant = 'atp'
                    energy_product = 'adp'
                elif compartment.id=='x':
                    energy_compartment = peroxisome
                    energy_reactant = 'atp'
                    energy_product = 'adp'
                else:
                    energy_compartment = er
                    energy_reactant = 'gtp'
                    energy_product = 'gdp'            

                # Adding participants to LHS
                # Include ATP/GTP hydrolysis during (co-translational and post-translational) translocation
                trans_reaction.participants.append(protein_model.species_coefficients.get_or_create(
                    coefficient=-1))
                trans_reaction.participants.append(model.species_types.get_one(id=energy_reactant).species.get_one(
                    compartment = energy_compartment).species_coefficients.get_or_create(coefficient=-1))
                trans_reaction.participants.append(model.species_types.get_one(id='h2o').species.get_one(
                    compartment = energy_compartment).species_coefficients.get_or_create(coefficient=-1))

                # Adding participants to RHS
                trans_reaction.participants.append(protein_model.species_type.species.get_one(
                    compartment=compartment).species_coefficients.get_or_create(coefficient=1))
                trans_reaction.participants.append(model.species_types.get_one(id=energy_product).species.get_one(
                    compartment = energy_compartment).species_coefficients.get_or_create(coefficient=1))
                trans_reaction.participants.append(model.species_types.get_one(id='pi').species.get_one(
                    compartment = energy_compartment).species_coefficients.get_or_create(coefficient=1))
                trans_reaction.participants.append(model.species_types.get_one(id='h').species.get_one(
                    compartment = energy_compartment).species_coefficients.get_or_create(coefficient=1))
                
                trans_rxn_no += 1

        print('{} reactions each for initiation and elongation and {} reactions for protein translocation '
            'have been generated'.format(init_el_rxn_no, trans_rxn_no))        

    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """                     
        cytoplasmic_ribosome = self.options.get('cytoplasmic_ribosome')
        mitochondrial_ribosome = self.options.get('mitochondrial_ribosome')
        cytoplasmic_initiation_factors = self.options.get('cytoplasmic_initiation_factors')
        mitochondrial_initiation_factors = options.get('mitochondrial_initiation_factors')
        cytoplasmic_elongation_factors = options.get('cytoplasmic_elongation_factors')
        mitochondrial_elongation_factors = options.get('mitochondrial_elongation_factors')
        cytoplasmic_chaperones = options.get('cytoplasmic_chaperones')
        mitochondrial_chaperones = options.get('mitochondrial_chaperones')
        er_chaperones = options.get('er_chaperones')
        
        cytosol = model.compartments.get_one(id='c')
        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')

        max_bool = model.parameters.get_or_create(
            id='max_bool_substance',
            type=None,
            value=1,
            units=unit_registry.parse_units('molecule'),
            comments='Boolean switch for determining if binding site is still available'
            )

        min_bool = model.parameters.get_or_create(
            id='min_bool_substance',
            type=None,
            value=0,
            units=unit_registry.parse_units('molecule'),
            comments='Boolean switch for determining if gene binding site is still available'
            )      

        # Generate response function for the tRNA(s) of each codon and for each amino acid
        trna_kb = cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType,
                                        type=wc_kb.eukaryote.TranscriptType.tRna)
        trna_grouping = {'c': {}, 'm': {}}
        for trna in trna_kb:
            anticodon_prop = i.properties.get_one(property='anticodon:amino_acid').get_value().split(':')
            codon = Bio.Seq.transcribe(Bio.Seq.reverse_complement(anticodon_prop[0]))
            if trna.species[0].compartment.id == 'm':
                if codon in trna_grouping['m']: 
                    trna_grouping['m'][codon]['trna'].append(trna.id)
                else:
                    trna_grouping['m'][codon] = {'trna': [trna.id], 'aa': anticodon_prop[1]}
            else:
                if codon in trna_grouping['c']: 
                    trna_grouping['c'][codon]['trna'].append(trna.id)
                else:
                    trna_grouping['c'][codon] = {'trna': [trna.id], 'aa': anticodon_prop[1]}
        
        trna_functions = {'c': {}, 'm': {}}
        for comp, all_trnas in trna_grouping.items():
            for codon, trnas in all_trnas.items():
                compartment = mitochondrion if comp=='m' else nucleus
                factor_exp, all_species, all_parameters, all_volumes, all_observables = self.gen_response_functions(
                    'translation', compartment, [trnas['trna']])

                objects = {
                    wc_lang.Species: all_species,
                    wc_lang.Parameter: all_parameters,
                    wc_lang.Observable: all_observables,
                    wc_lang.Function: all_volumes,            
                    }

                trna_expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                assert error is None, str(error)

                trna_functions[comp][codon] = {
                    'function': model.functions.create(     
                            id='trna_function_{}_{}'.format(codon, compartment.id),               
                            name='tRNA response function for codon {} in {}'.format(codon, compartment.name),
                            expression=trna_expression,
                            units=unit_registry.parse_units(''),
                            ),
                    'aa': trnas['aa'],
                    'objects':objects,
                    }

        aa_functions = {'c': {}, 'm': {}}
        for aa_id in amino_acid_id_conversion.values():
            for compartment in [cytosol, mitochondrion]:
                factor_exp, all_species, all_parameters, all_volumes, all_observables = self.gen_response_functions(
                    'translation', compartment, [aa_id])

                objects = {
                    wc_lang.Species: all_species,
                    wc_lang.Parameter: all_parameters,
                    wc_lang.Observable: all_observables,
                    wc_lang.Function: all_volumes,            
                    }                                
                
                aa_expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                assert error is None, str(error)

                aa_functions[compartment.id] = {
                    aa_id: {
                        'function': model.functions.create(     
                            id='aminoacid_function_{}_{}'.format(aa_id, compartment.id),               
                            name='response function for amino acid {} in {}'.format(aa_id, compartment.name),
                            expression=aa_expression,
                            units=unit_registry.parse_units(''),
                            ),
                        'objects': objects},
                    }            
        
        rate_law_no = 0  
        mrna_kbs = [i for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType) \
            if i.type==wc_kb.eukaryote.TranscriptType.mRna]      
        for mrna_kb in mrna_kbs:
            
            mrna_kb_compartment_id = mrna_kb.species[0].compartment.id
            if mrna_kb_compartment_id == 'n':
                ribosome_complex = model.species_types.get_one(
                    name=cytoplasmic_ribosome).species.get_one(compartment=cytosol)
                initiation_factors = cytoplasmic_initiation_factors
                elongation_factors = cytoplasmic_elongation_factors
                translation_compartment = cytosol 
            else:
                ribosome_complex = model.species_types.get_one(
                    name=mitochondrial_ribosome).species.get_one(compartment=mitochondrion)
                initiation_factors = mitochondrial_initiation_factors
                elongation_factors = mitochondrial_elongation_factors
                translation_compartment = mitochondrion

            # Generate rate law for initiation
            init_reaction = model.reactions.get_one(id='translation_initiation_' + mrna_kb.id)
            factor_exp, all_species, all_parameters, all_volumes, all_observables = self.gen_response_functions(
                init_reaction.id, translation_compartment, initiation_factors)
            
            specific_binding_constant = model.parameters.create(
                id='{}_ribosome_binding_constant'.format(mrna_kb.id),
                type=None,
                units=unit_registry.parse_units('molecule^-2 s^-1'),
                )
            
            all_species[ribosome_complex.id] = ribosome_complex
            all_species[self._allowable_queue_len[mrna_kb.id][0].id]= self._allowable_queue_len[mrna_kb.id][0]

            all_parameters[max_bool.id] = max_bool
            all_parameters[min_bool.id] = min_bool
            all_parameters[specific_binding_constant.id] = specific_binding_constant
            
            expression = '{} * {} * max(min({} , {}) , {}){}'.format(
                specific_binding_constant.id,
                ribosome_complex.id,
                self._allowable_queue_len[mrna_kb.id][0].id,
                max_bool.id,
                min_bool.id,
                (' * {} * 2**{}'.format(' * '.join(factor_exp), len(factor_exp))) if factor_exp else '',                
                )

            init_rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, {
                wc_lang.Species: all_species,
                wc_lang.Parameter: all_parameters,
                wc_lang.Observable: all_observables,
                wc_lang.Function: all_volumes,            
                })
            assert error is None, str(error)

            init_rate_law = model.rate_laws.create(
                direction=wc_lang.RateLawDirection.forward,
                type=None,
                expression=init_rate_law_expression,
                reaction=init_reaction,
                units=unit_registry.parse_units('s^-1'),                
                )
            init_rate_law.id = init_rate_law.gen_id()
            
            # Generate rate law for elongation and termination
            elongation_reaction = model.reactions.get_one(id='translation_elongation_' + mrna_kb.id)
            expression_terms, all_species, all_parameters, all_volumes, all_observables = self.gen_response_functions(
                elongation_reaction.id, translation_compartment, elongation_factors + [['gtp'], ['atp']])
            
            ribo_bound_species = model.species_types.get_one(id='ribo_bound_{}'.format(mrna_kb.id)).species[0]
            all_species[ribo_bound_species.id] = ribo_bound_species

            k_cat_elongation = model.parameters.create(
                id='k_cat_{}'.format(elongation_reaction.id),
                type=wc_ontology['WC:k_cat'],
                units=unit_registry.parse_units('molecule^-1 s^-1'),
                )
            all_parameters[k_cat_elongation.id] = k_cat_elongation

            objects = {
                wc_lang.Species: all_species,
                wc_lang.Parameter: all_parameters,
                wc_lang.Observable: all_observables,
                wc_lang.Function: all_volumes,            
                }
            
            mrna_seq = mrna_kb.get_seq()
            all_codons = set([mrna_seq[i * 3:(i + 1) * 3] for i in range((len(mrna_seq) + 3 - 1) // 3 )])            
            for i in all_codons:                
                if len(i)==3 and i in trna_functions[translation_compartment.id]:

                    codon_info = trna_functions[translation_compartment.id][i]
                    expression_terms.append(codon_info['function'].expression.expression)
                    
                    for cl, dictionary in objects.items():
                        dictionary.update(codon_info['objects'][cl])
                    
                    if aa_functions[translation_compartment.id][codon_info['aa']][
                        'function'].expression.expression not in expression_terms:
                        
                        expression_terms.append(aa_functions[translation_compartment.id][
                            codon_info['aa']]['function'].expression.expression)

                        for cl, dictionary in objects.items():
                            dictionary.update(aa_functions[translation_compartment.id][
                                codon_info['aa']]['objects'][cl])
            
            expression = '{} * {} * {} * 2**{}'.format(
                k_cat_elongation.id,
                ribo_bound_species.id,
                ' * '.join(expression_terms),
                len(expression_terms),
                )

            el_rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, objects)
            assert error is None, str(error)
            
            el_rate_law = model.rate_laws.create(
                direction=wc_lang.RateLawDirection.forward,
                type=None,
                expression=el_rate_law_expression,
                reaction=elongation_reaction,
                )
            el_rate_law.id = el_rate_law.gen_id()
            
            rate_law_no += 1
        
        # Generate rate law for translocation
        trans_rxn_no = 0    
        for reaction in self.submodel.reactions:
            if 'translocation' in reaction.id:

                translation_compartment = model.compartments.get_one(id=reaction.id.split('_')[2])                
                target_compartment_id = '_'.join(reaction.id.split('_')[4:])
                
                if target_compartment_id=='n'
                    energy_compartment = nucleus
                    energy_reactant = 'gtp'       
                    chaperones = []
                elif target_compartment_id=='m':
                    energy_compartment = mitochondrion
                    energy_reactant = 'atp'
                    chaperones = mitochondrial_chaperones                    
                elif target_compartment_id=='x':
                    energy_compartment = peroxisome
                    energy_reactant = 'atp'
                    chaperones = []
                else:
                    energy_compartment = er
                    energy_reactant = 'gtp'
                    chaperones = er_chaperones

                factor_exp, all_species, all_parameters, all_volumes, all_observables = self.gen_response_functions(
                    reaction.id, energy_compartment, chaperones + [['gtp'], ['atp']])

                k_cat_translocation = model.parameters.create(
                    id='k_cat_{}'.format(reaction.id),
                    type=wc_ontology['WC:k_cat'],
                    units=unit_registry.parse_units('molecule^-1 s^-1'),
                    )
                all_parameters[k_cat_translocation.id] = k_cat_translocation

                protein_species = model.species_types.get_one(
                    id=reaction.id.split('_')[1]).species.get_one(compartment=translation_compartment)
                all_species[protein_species.id] = protein_species

                volume = translation_compartment.init_density.function_expressions[0].function
                all_volumes[volume.id] = volume

                expression = '{} * {} * {} * 2**{}'.format(
                    k_cat_translocation.id,
                    protein_species.id,
                    ' * '.join(factor_exp),
                    len(factor_exp),
                    )

                trans_rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, {
                    wc_lang.Species: all_species,
                    wc_lang.Parameter: all_parameters,
                    wc_lang.Observable: all_observables,
                    wc_lang.Function: all_volumes,            
                    })
                assert error is None, str(error)
                
                trans_rate_law = model.rate_laws.create(
                    direction=wc_lang.RateLawDirection.forward,
                    type=None,
                    expression=trans_rate_law_expression,
                    reaction=reaction,
                    )
                trans_rate_law.id = trans_rate_law.gen_id()
                    
                trans_rxn_no += 1
                    
        print('{} rate laws for initiation reactions, {} rate laws for elongation '
            'reactions and {} rate laws for translocation reactions have been generated'.format(
            rate_law_no, rate_law_no, trans_rxn_no))            
                        
    def calibrate_submodel(self):
        """ Calibrate the submodel using data in the KB """
        
        model = self.model        
        cell = self.knowledge_base.cell

        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')
        cytosol = model.compartments.get_one(id='c')

        init_compartment_volumes = {
            nucleus.id: nucleus.init_volume.mean * nucleus.init_density.value,
            mitochondrion.id: mitochondrion.init_volume.mean * mitochondrion.init_density.value,
            cytosol.id: cytosol.init_volume.mean * cytosol.init_density.value,
            }

        beta = self.options.get('beta')
        polysome_fraction = self.options['polysome_fraction']
        cytoplasmic_ribosome = self.options.get('cytoplasmic_ribosome')
        mitochondrial_ribosome = self.options.get('mitochondrial_ribosome')

        cytoplasmic_ribosome_species = model.species_types.get_one(
            name=cytoplasmic_ribosome).species.get_one(compartment=cytosol)
        mitochondrial_ribosome_species = model.species_types.get_one(
            name=mitochondrial_ribosome).species.get_one(compartment=mitochondrion)

        Avogadro = self.model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))

        mean_doubling_time = model.parameters.get_one(id='mean_doubling_time').value       
        
        mrna_kbs = [i for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType) \
            if i.type==wc_kb.eukaryote.TranscriptType.mRna]
        
        # Determine initial concentrations of ribosome bound sites and update the concentrations of free ribosomes
        cytoplasmic_bound_ribosomes = 0
        mitochondrial_bound_ribosomes = 0        
        for mrna_kb in mrna_kbs:            

            mrna_kb_compartment_id = mrna_kb.species[0].compartment.id
            
            if mrna_kb_compartment_id == 'n':   

                ribo_bound_conc = polysome_fraction[mrna_kb.id] * \
                    cytoplasmic_ribosome_species.distribution_init_concentration.mean                
                cytoplasmic_bound_ribosomes += ribo_bound_conc
            
            else:

                ribo_bound_conc = polysome_fraction[mrna_kb.id] * \
                    mitochondrial_ribosome_species.distribution_init_concentration.mean
                mitochondrial_bound_ribosomes += ribo_bound_conc

            ribo_bound_species = model.species_types.get_one(id='ribo_bound_{}'.format(
                mrna_kb.id)).species[0]             
            ribo_bound_species.distribution_init_concentration.mean = ribo_bound_conc
        
        cytoplasmic_ribosome_species.distribution_init_concentration.mean -= cytoplasmic_bound_ribosomes    
        mitochondrial_ribosome_species.distribution_init_concentration.mean -= mitochondrial_bound_ribosomes
             
        # Calibrate initiation and elongation reactions
        undetermined_el_km = []
        determined_el_km = []
        undetermined_trans_km = []
        determined_trans_km = []              
        for mrna_kb in mrna_kbs:

            mrna_kb_compartment_id = mrna_kb.species[0].compartment.id
            ribosome_complex = cytoplasmic_ribosome_species if mrna_kb_compartment_id == 'n' else mitochondrial_ribosome_species

            protein_model = model.species_types.get_one(id=mrna_kb.protein.id)
            
            complex_model_stoic = {model.species_types.get_one(id=i.id):j.coefficient for i in cell.species_types.get(
                __type=wc_kb.core.ComplexSpeciesType) for j in i.subunits if j.species_type==mrna_kb.protein}

            total_concentration = sum([i.distribution_init_concentration.mean for i in protein_model.species]) + \
                sum([i.distribution_init_concentration.mean*v for k,v in complex_model_stoic.items() for i in k.species])
            half_life = mrna_kb.protein.properties.get_one(property='half-life').get_value()

            average_rate = utils.calc_avg_syn_rate(
                total_concentration, half_life, mean_doubling_time)

            # Calibrate initiation reaction
            init_reaction = model.reactions.get_one(id='translation_initiation_' + mrna_kb.id)
            
            init_species_counts[ribosome_complex.id] = ribosome_complex.distribution_init_concentration.mean
            init_species_counts[self._allowable_queue_len[mrna_kb.id][0].id] = self._allowable_queue_len[mrna_kb.id][1]

            model_kcat = model.parameters.get_one(id='{}_ribosome_binding_constant'.format(mrna_kb.id))

            if average_rate:            
                model_kcat.value = 1.
                model_kcat.value = average_rate / init_reaction.rate_laws[0].expression._parsed_expression.eval({
                    wc_lang.Species: init_species_counts,
                    wc_lang.Compartment: init_compartment_volumes,
                })
            else:          
                model_kcat.value = 0.

            # Calibrate elongation reaction
            el_reaction = model.reactions.get_one(id='translation_elongation_' + mrna_kb.id)
            
            el_species_counts = {}
            el_observable_values = {}            

            for species in el_reaction.rate_laws[0].expression.species:
                el_species_counts[species.id] = species.distribution_init_concentration.mean
                if model.parameters.get(id='K_m_{}_{}'.format(el_reaction.id, species.species_type.id)):
                    model_Km = model.parameters.get_one(
                        id='K_m_{}_{}'.format(el_reaction.id, species.species_type.id))
                    if species.distribution_init_concentration.mean:
                        model_Km.value = beta * species.distribution_init_concentration.mean \
                            / Avogadro.value / species.compartment.init_volume.mean
                        model_Km.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                            beta, species.species_type.id, species.compartment.name)
                        determined_el_km.append(model_Km.value)
                    else:
                        undetermined_el_km.append(model_Km)

            for obs in el_reaction.rate_laws[0].expression.observables:
                
                obs_species_counts = {}    
                for species in obs.species:
                    obs_species_counts[species.id] = species.distribution_init_concentration.mean
                    el_species_counts[species.id] = species.distribution_init_concentration.mean
                
                obs_value = obs.expression._parsed_expression.eval({
                    wc_lang.Species: obs_species_counts,
                    wc_lang.Compartment: init_compartment_volumes})
                el_observable_values[obs.id] = obs_value

                if model.parameters.get(id='K_m_{}_{}'.format(el_reaction.id, obs.id)):
                    model_Km = model.parameters.get_one(
                        id='K_m_{}_{}'.format(el_reaction.id, obs.id))
                    model_Km.value = beta * obs_value / Avogadro.value \
                            / species.compartment.init_volume.mean # All species in the obs should be in the same compartment
                    model_Km.comments = 'The value was assumed to be {} times the value of {}'.format(
                        beta, obs.id)                

            model_kcat = model.parameters.get_one(id='k_cat_{}'.format(el_reaction.id))

            if average_rate:            
                model_kcat.value = 1.
                model_kcat.value = average_rate / el_reaction.rate_laws[0].expression._parsed_expression.eval({
                    wc_lang.Species: el_species_counts,
                    wc_lang.Compartment: init_compartment_volumes,
                    wc_lang.Observable: el_observable_values,
                })
            else:          
                model_kcat.value = 0.

            # Calibrate translocation reaction
            conc_per_comp = {}
            for protein in protein_model.species:
                conc_per_comp[protein.compartment] = protein.distribution_init_concentration.mean
            for cplx_st, stoic in complex_model_stoic.items():
                for cplx_species in cplx_st.species:
                    if cplx_species.compartment in conc_per_comp:
                        conc_per_comp[cplx_species.compartment] += stoic * \
                            cplx_species.distribution_init_concentration.mean
                    else:        
                        conc_per_comp[cplx_species.compartment] = stoic * \
                            cplx_species.distribution_init_concentration.mean
            
            translation_compartment = cytosol if mrna_kb_compartment_id == 'n' else mitochondrion

            for compartment, trans_reaction in self._translocation_reactions[mrna_kb]:
                trans_species_counts = {}
                trans_observable_values = {}            

                for species in trans_reaction.rate_laws[0].expression.species:
                    trans_species_counts[species.id] = species.distribution_init_concentration.mean
                    if model.parameters.get(id='K_m_{}_{}'.format(trans_reaction.id, species.species_type.id)):
                        model_Km = model.parameters.get_one(
                            id='K_m_{}_{}'.format(trans_reaction.id, species.species_type.id))
                        if species.distribution_init_concentration.mean:
                            model_Km.value = beta * species.distribution_init_concentration.mean \
                                / Avogadro.value / species.compartment.init_volume.mean
                            model_Km.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                                beta, species.species_type.id, species.compartment.name)
                            determined_trans_km.append(model_Km.value)
                        else:
                            undetermined_trans_km.append(model_Km)

                for obs in trans_reaction.rate_laws[0].expression.observables:
                    
                    obs_species_counts = {}    
                    for species in obs.species:
                        obs_species_counts[species.id] = species.distribution_init_concentration.mean
                        trans_species_counts[species.id] = species.distribution_init_concentration.mean
                    
                    obs_value = obs.expression._parsed_expression.eval({
                        wc_lang.Species: obs_species_counts,
                        wc_lang.Compartment: init_compartment_volumes})
                    trans_observable_values[obs.id] = obs_value

                    if model.parameters.get(id='K_m_{}_{}'.format(trans_reaction.id, obs.id)):
                        model_Km = model.parameters.get_one(
                            id='K_m_{}_{}'.format(trans_reaction.id, obs.id))
                        model_Km.value = beta * obs_value / Avogadro.value \
                                / species.compartment.init_volume.mean # All species in the obs should be in the same compartment
                        model_Km.comments = 'The value was assumed to be {} times the value of {}'.format(
                            beta, obs.id)                

                model_kcat = model.parameters.get_one(id='k_cat_{}'.format(trans_reaction.id))

                if average_rate:            
                    model_kcat.value = 1.
                    model_kcat.value = conc_per_comp[compartment] / conc_per_comp[translation_compartment] * \
                        average_rate / trans_reaction.rate_laws[0].expression._parsed_expression.eval({
                        wc_lang.Species: trans_species_counts,
                        wc_lang.Compartment: init_compartment_volumes,
                        wc_lang.Observable: trans_observable_values,
                        })
                else:          
                    model_kcat.value = 0.                
        
        median_km = numpy.median(determined_el_km)
        for model_Km in undetermined_el_km:
            model_Km.value = median_km
            model_Km.comments = 'Set to the median value because transcript concentration was zero'

        median_km = numpy.median(determined_trans_km)
        for model_Km in undetermined_trans_km:
            model_Km.value = median_km
            model_Km.comments = 'Set to the median value because transcript concentration was zero'    

    def gen_response_functions(self, reaction_id, compartment, reaction_factors):
        """ Generate a list of response function expression string for each factor or 
            group of factors (F) in the form of:
                       
                        F/(Km + F)

        Args:
            reaction_id (:obj:`str`): identifier of reaction whose rate law will use the function expressions
            compartment (:obj:`wc_lang.Compartment`): compartment where the reaction occurs
            reaction_factors (:obj:`list` of `list`): list of lists of the name of
                (initiation/elongation/translocation) factors, grouped based on similar functions or classes
            
        Returns:
            :obj:`list`: list of response function expression string for each factor/group of factors
            :obj:`dict`: IDs of species (keys) and their species objects (values)
            :obj:`dict`: IDs of parameters (keys) and their parameter objects (values)
            :obj:`dict`: IDs of volume density functions (keys) and their function objects (values)
            :obj:`dict`: IDs of observables (keys) and their observable objects (values)
        """
        model = self.model

        all_species = {}
        all_parameters = {}
        all_volumes = {}
        all_observables = {}

        avogadro = model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))
        all_parameters[avogadro.id] = avogadro

        volume = compartment.init_density.function_expressions[0].function
        all_volumes[volume.id] = volume

        factor_exp = []
        n = 0
        for factors in reaction_factors:
            
            if len(factors) == 1:
                
                factor_species = model.species_types.get_one(
                    name=factors[0]).species.get_one(compartment=compartment)
                all_species[factor_species.gen_id()] = factor_species

                model_k_m = model.parameters.create(
                    id='K_m_{}_{}'.format(reaction_id, factor_species.species_type.id),
                    type=wc_ontology['WC:K_m'],
                    units=unit_registry.parse_units('M'))
                all_parameters[model_k_m.id] = model_k_m                    

                factor_exp.append('({} / ({} + {} * {} * {}))'.format(
                    factor_species.gen_id(),
                    factor_species.gen_id(),
                    model_k_m.id, 
                    avogadro.id,
                    volume.id))

            else:
                
                obs_exp = []                    
                for factor in factors:                        
                    factor_species = model.species_types.get_one(
                        name=factors).species.get_one(compartment=compartment)
                    all_species[factor_species.gen_id()] = factor_species
                    obs_exp.append(factor_species.gen_id())
                
                n += 1
                
                observable_exp, error = wc_lang.ObservableExpression.deserialize(
                ' + '.join(obs_exp),
                {wc_lang.Species: all_species})            
                assert error is None, str(error)                
                
                factor_observable = model.observables.create(
                    id='{}_factors_{}_{}'.format(reaction_id, compartment.id, n), 
                    name='factor for {} in {} {}'.format(reaction_id, compartment.name, n), 
                    units=unit_registry.parse_units('molecule'), 
                    expression=observable_exp)
                all_observables[factor_observable.id] = factor_observable

                model_k_m = model.parameters.create(
                    id='K_m_{}_{}'.format(reaction_id, factor_observable.id),
                    type=wc_ontology['WC:K_m'],
                    units=unit_registry.parse_units('M'))
                all_parameters[model_k_m.id] = model_k_m

                factor_exp.append('({} / ({} + {} * {} * {}))'.format(
                    factor_observable.id,
                    factor_observable.id,
                    model_k_m.id, 
                    avogadro.id,
                    volume.id))

        return factor_exp, all_species, all_parameters, all_volumes, all_observables            
