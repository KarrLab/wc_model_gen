""" Generator for transcription submodels for eukaryotes

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-07
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_onto import onto
from wc_utils.util.chem import EmpiricalFormula
from wc_utils.util.units import unit_registry
import wc_model_gen.global_vars as gvar
import wc_model_gen.utils as utils
import collections
import math
import numpy
import scipy.constants
import wc_kb
import wc_lang
import wc_model_gen


class TranscriptionSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for transcription submodel 

        Options:
        * rna_pol_pair (:obj:`dict`): a dictionary of RNA id as key and
            the name of RNA polymerase complex that transcribes the RNA as value, e.g.
            rna_pol_pair = {
                'rRNA45S': 'DNA-directed RNA Polymerase I complex', 
                'mRNA': 'DNA-directed RNA Polymerase II complex',
                'sRNA': 'DNA-directed RNA Polymerase II complex', 
                'tRNA': 'DNA-directed RNA Polymerase III complex',
                'rRNA5S': 'DNA-directed RNA Polymerase III complex'
                }
        * init_factors (:obj:`dict`, optional): a dictionary of generic init factor name as key and 
            list of lists of the id or name of initiation factors, grouped based on similar functions or classes,
            e.g. {'pol1_init_factors': [['factor1_variant1', 'factor1_variant2'], ['factor2']]} where the keys must start with
            the substring 'pol1_', 'pol2_', 'pol3_', and 'polm_' if factors for the polymerase exists, 
            the default is an empty dictionary
        * elongation_termination_factors (:obj:`dict`, optional): a dictionary of generic elongation and 
            termination factor name as key and list of lists of the id or name of elongation and termination factors, 
            grouped based on similar functions or classes, 
            e.g. {'pol1_elongation_termination_factors': [['factor1_variant1', 'factor1_variant2'], ['factor2']]},
            where the keys must start with the substring 'pol1_', 'pol2_', 'pol3_', and 'polm_' if factors for the polymerase exists,
            the default is an empty dictionary
        * elongation_negative_factors (:obj:`dict`, optional): a dictionary of generic elongation negative
            factor name as key and list of lists of the id or name of elongation negative factors, 
            grouped based on similar functions or classes, 
            e.g. {'pol2_elongation_negative_factors': [['factor1_variant1', 'factor1_variant2'], ['factor2']]},
            where the keys must start with the substring 'pol1_', 'pol2_', 'pol3_', and 'polm_' if factors for the polymerase exists,
            the default is an empty dictionary       
        * rna_init_factors (:obj:`dict`, optional): a dictionary of RNA id as key and the generic init factor
            name (the key in init_factors option) as value, the default is an empty dictionary   
        * rna_elongation_termination_factors (:obj:`dict`, optional): a dictionary of RNA id as key and the 
            generic elongation and termination factor name (the key in elongation_termination_factors option) 
            as value, the default is an empty dictionary
        * rna_elongation_negative_factors (:obj:`dict`, optional): a dictionary of RNA id as key and the 
            generic elongation negatic factor name (the key in elongation_termination_factors option) as value,
            the default is an empty dictionary        
        * beta (:obj:`float`, optional): ratio of Michaelis-Menten constant to substrate 
            concentration (Km/[S]) for use when estimating Km values, the default value is 1
        * beta_activator (:obj:`float`, optional): ratio of effective equilibrium 
            dissociation constant of a transcription factor (activator) to the transcription 
            factor concentration (Ka/[TF]) for use when estimating Ka values, the default value is 1
        * beta_repressor (:obj:`float`, optional): ratio of effective equilibrium 
            dissociation constant of a transcription factor (repressor) to the transcription 
            factor concentration (Kr/[TF]) for use when estimating Kr values, the default value is 1
        * activator_effect (:obj:`float`, optional): interaction effect between an activator 
            and RNA polymerase, which must take the value of 1 and higher, the default value is 1.2    
        * polr_occupancy_width (:obj:`int`, optional): number of base-pairs on the DNA occupied 
            by each bound RNA polymerase, , the default value is 80
        * ribosome_occupancy_width (:obj:`int`, optional): number of base-pairs on the mRNA occupied 
            by each bound ribosome, the default value is 27 (9 codons)                              
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        if 'rna_pol_pair' not in options:
            raise ValueError('The dictionary rna_pol_pair has not been provided')
        else:    
            rna_pol_pair = options['rna_pol_pair']

        init_factors = options.get('init_factors', {})
        options['init_factors'] = init_factors
        if init_factors:
            if not any(j in i for i in init_factors.keys() for j in ['pol1_', 'pol2_', 'pol3_', 'polm_']):
                raise ValueError(
                    'The keys in init_factors must start with the substrings "pol1_", "pol2_", "pol3_" and/or "polm_"')

        elongation_termination_factors = options.get('elongation_termination_factors', {})
        options['elongation_termination_factors'] = elongation_termination_factors
        if elongation_termination_factors:
            if not any(j in i for i in elongation_termination_factors.keys() for j in ['pol1_', 'pol2_', 'pol3_', 'polm_']):
                raise ValueError(
                    'The keys in elongation_termination_factors must start with the substrings "pol1_", "pol2_", "pol3_" and/or "polm_"')

        elongation_negative_factors = options.get('elongation_negative_factors', {})
        options['elongation_negative_factors'] = elongation_negative_factors
        if elongation_negative_factors:
            if not any(j in i for i in elongation_negative_factors.keys() for j in ['pol1_', 'pol2_', 'pol3_', 'polm_']):
                raise ValueError(
                    'The keys in elongation_negative_factors must start with the substrings "pol1_", "pol2_", "pol3_" and/or "polm_"')

        rna_init_factors = options.get('rna_init_factors', {})
        options['rna_init_factors'] = rna_init_factors
        if rna_init_factors:
            for i in set(rna_init_factors.values()):
                if i not in init_factors:
                    raise ValueError('{} is not a key in init_factors'.format(i))

        rna_elongation_termination_factors = options.get('rna_elongation_termination_factors', {})
        options['rna_elongation_termination_factors'] = rna_elongation_termination_factors
        if rna_elongation_termination_factors:
            for i in set(rna_elongation_termination_factors.values()):
                if i not in elongation_termination_factors:
                    raise ValueError('{} is not a key in elongation_termination_factors'.format(i))

        rna_elongation_negative_factors = options.get('rna_elongation_negative_factors', {})
        options['rna_elongation_negative_factors'] = rna_elongation_negative_factors
        if rna_elongation_negative_factors:
            for i in set(rna_elongation_negative_factors.values()):
                if i and i not in elongation_negative_factors:
                    raise ValueError('{} is not a key in elongation_negative_factors'.format(i))    

        beta = options.get('beta', 1.)
        options['beta'] = beta

        beta_activator = options.get('beta_activator', 1.)
        options['beta_activator'] = beta_activator

        beta_repressor = options.get('beta_repressor', 1.)
        options['beta_repressor'] = beta_repressor

        activator_effect = options.get('activator_effect', 1.2)
        options['activator_effect'] = activator_effect

        polr_occupancy_width = options.get('polr_occupancy_width', 80)
        options['polr_occupancy_width'] = polr_occupancy_width

        ribosome_occupancy_width = options.get('ribosome_occupancy_width', 27)
        options['ribosome_occupancy_width'] = ribosome_occupancy_width

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')
        cytosol = model.compartments.get_one(id='c')

        polr_occupancy_width = self.options.get('polr_occupancy_width')
        ribosome_occupancy_width = self.options['ribosome_occupancy_width']        

        # Get species involved in reaction
        metabolic_participants = ['atp', 'ctp', 'gtp', 'utp', 'ppi', 
            'amp', 'cmp', 'gmp', 'ump', 'h2o', 'h', 'adp', 'pi']
        metabolites = {}
        for met in metabolic_participants:
            met_species_type = model.species_types.get_one(id=met)
            metabolites[met] = {
                'n': met_species_type.species.get_or_create(compartment=nucleus, model=model),
                'm': met_species_type.species.get_or_create(compartment=mitochondrion, model=model)
                }

        ref_polr_width = wc_lang.Reference(
            model=model,
            title='Structure and mechanism of the RNA Polymerase II transcription machinery',
            author='Steven Hahn',        
            year=2004,
            type=onto['WC:article'],
            publication='Nature Structural & Molecular Biology',        
            volume='11',
            issue='5',
            pages='394-403'
            )
        ref_polr_width.id = 'ref_'+str(len(model.references))

        ref_polr_distribution = wc_lang.Reference(
            model=model,
            title='In vivo dynamics of RNA polymerase II transcription',
            author='Xavier Darzacq, Yaron Shav-Tal, Valeria de Turris, Yehuda Brody, '
                'Shailesh M Shenoy, Robert D Phair, Robert H Singer',        
            year=2007,
            type=onto['WC:article'],
            publication='Nature Structural & Molecular Biology',        
            volume='14',
            pages='796-806'
            )
        ref_polr_distribution.id = 'ref_'+str(len(model.references))

        ref_ribo_width = wc_lang.Reference(
            model=model,
            title='Genome-wide analysis in vivo of translation with nucleotide resolution using ribosome profiling',
            author='Nicholas T Ingolia, Sina Ghaemmaghami, John R. S. Newman, Jonathan S. Weissman',        
            year=2009,
            type=onto['WC:article'],
            publication='Science',        
            volume='324',
            issue='5924',
            pages='218-223'
            )
        ref_ribo_width.id = 'ref_'+str(len(model.references))  

        print('Start generating transcription submodel...')                  
        
        # Create for each RNA polymerase a reaction of binding to non-specific site        
        nuclear_genome_length = 0
        mitochondrial_genome_length = 0
        for chromosome in cell.species_types.get(__type=wc_kb.core.DnaSpeciesType):
            if 'M' in chromosome.id:
                mitochondrial_genome_length += len(chromosome.get_seq())
            else:
                nuclear_genome_length += len(chromosome.get_seq())
        self._mitochondrial_max_binding_sites = math.floor(
            mitochondrial_genome_length/polr_occupancy_width)
        self._nuclear_max_binding_sites = math.floor(
            nuclear_genome_length/polr_occupancy_width)

        self._total_polr = {}
        self._gene_bound_polr = {}
        rna_pol_pair = self.options.get('rna_pol_pair')
        for polr in set(rna_pol_pair.values()):

            self._gene_bound_polr[polr] = []
            
            if 'mito' in polr:
                transcription_compartment = mitochondrion
                genome_sites = self._mitochondrial_max_binding_sites
            else:
                transcription_compartment = nucleus
                genome_sites = self._nuclear_max_binding_sites
            
            polr_complex = model.species_types.get_one(name=polr)
            polr_complex_species = model.species.get_one(
                species_type=polr_complex, compartment=transcription_compartment)
            conc_free_polr = model.distribution_init_concentrations.get_one(
                species=polr_complex_species)
            self._total_polr[polr] = conc_free_polr.mean
            conc_free_polr.mean = math.floor(0.75*conc_free_polr.mean)
            conc_free_polr.comments = 'The free pool is estimated to be three quarters of the total concentration'
            conc_free_polr.references.append(ref_polr_distribution)
            
            polr_non_specific_binding_site_st = model.species_types.get_or_create(
                id='polr_non_specific_binding_site',
                name='non-specific binding site of RNA polymerase',
                type=onto['WC:pseudo_species'],
                )
            polr_non_specific_binding_site_st.structure = wc_lang.ChemicalStructure(
                    empirical_formula = EmpiricalFormula(),
                    molecular_weight = 0.,
                    charge = 0)
            polr_non_specific_binding_site_species = model.species.get_or_create(
                species_type=polr_non_specific_binding_site_st, compartment=transcription_compartment)
            polr_non_specific_binding_site_species.id = polr_non_specific_binding_site_species.gen_id()

            conc_model = model.distribution_init_concentrations.get_or_create(
                species=polr_non_specific_binding_site_species,
                mean=genome_sites,
                units=unit_registry.parse_units('molecule'),
                comments='Set to genome length divided by {} bp to allow '
                    'queueing of RNA polymerase during transcription'.format(polr_occupancy_width),
                references=[ref_polr_width],
                )
            conc_model.id = conc_model.gen_id()

            polr_bound_non_specific_species_type = model.species_types.get_or_create(
                id='{}_bound_non_specific_site'.format(polr_complex.id),
                name='{}-bound non-specific site'.format(polr_complex.id),
                type=onto['WC:pseudo_species'],                
                )
            polr_bound_non_specific_species_type.structure = wc_lang.ChemicalStructure(
                    empirical_formula = polr_complex.structure.empirical_formula,
                    molecular_weight = polr_complex.structure.molecular_weight,
                    charge = polr_complex.structure.charge)
            polr_bound_non_specific_species = model.species.get_or_create(
                species_type=polr_bound_non_specific_species_type, compartment=transcription_compartment)
            polr_bound_non_specific_species.id = polr_bound_non_specific_species.gen_id()

            conc_model = model.distribution_init_concentrations.get_or_create(
                species=polr_bound_non_specific_species,
                mean=math.floor(self._total_polr[polr]*0.2475),
                units=unit_registry.parse_units('molecule'),
                comments='Approximately 24.75 percent of RNA polymerase is bound to non-specific site',
                references=[ref_polr_distribution])
            conc_model.id = conc_model.gen_id()

            ns_binding_reaction = model.reactions.create(
                submodel=self.submodel, id='non_specific_binding_{}'.format(polr_complex.id),
                name='non-specific binding of {} in {}'.format(polr, transcription_compartment.name),
                reversible=False)
            
            ns_binding_reaction.participants.append(
                polr_complex_species.species_coefficients.get_or_create(
                coefficient=-1))
            ns_binding_reaction.participants.append(
                polr_non_specific_binding_site_species.species_coefficients.get_or_create(
                coefficient=-1))
            ns_binding_reaction.participants.append(
                polr_bound_non_specific_species.species_coefficients.get_or_create(
                coefficient=1))
        
        # Create initiation and elongation reactions for each RNA
        init_el_rxn_no = 0
        rna_kbs = cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType)
        self._initiation_polr_species = {}
        self._elongation_modifier = {}
        self._allowable_queue_len = {}        
        for rna_kb in rna_kbs:
            
            transcription_compartment = nucleus if rna_kb.species[0].compartment.id == 'c' else mitochondrion
            translation_compartment = cytosol if rna_kb.species[0].compartment.id == 'c' else mitochondrion
            
            # Create initiation reaction
            polr_complex = model.species_types.get_one(name=rna_pol_pair[rna_kb.id])
            polr_complex_species = model.species.get_one(
                species_type=polr_complex, compartment=transcription_compartment)
            self._initiation_polr_species[rna_kb.id] = polr_complex_species
            
            polr_bound_non_specific_species_type = model.species_types.get_one(
                id='{}_bound_non_specific_site'.format(polr_complex.id))
            polr_bound_non_specific_species = model.species.get_one(
                species_type=polr_bound_non_specific_species_type, compartment=transcription_compartment)
            
            polr_non_specific_binding_site_st = model.species_types.get_one(
                id='polr_non_specific_binding_site')
            polr_non_specific_binding_site_species = model.species.get_one(
                species_type=polr_non_specific_binding_site_st, compartment=transcription_compartment)
            
            gene = rna_kb.gene
            polr_binding_site_st = model.species_types.get_or_create(
                id='{}_binding_site'.format(gene.id),
                name='binding site of {}'.format(gene.name),
                type=onto['WC:pseudo_species'],
                )
            polr_binding_site_st.structure = wc_lang.ChemicalStructure(
                    empirical_formula = EmpiricalFormula(),
                    molecular_weight = 0.,
                    charge = 0)
            polr_binding_site_species = model.species.get_or_create(
                species_type=polr_binding_site_st, compartment=transcription_compartment)
            polr_binding_site_species.id = polr_binding_site_species.gen_id()

            gene_seq = gene.get_seq() 
            conc_model = model.distribution_init_concentrations.create(
                species=polr_binding_site_species,
                mean=math.floor(len(gene_seq)/polr_occupancy_width) + 1,
                units=unit_registry.parse_units('molecule'),
                comments='Set to gene length divided by {} bp to allow '
                    'queueing of RNA polymerase during transcription'.format(polr_occupancy_width),
                references=[ref_polr_width]    
                )
            conc_model.id = conc_model.gen_id()
            self._allowable_queue_len[rna_kb.id] = (polr_binding_site_species, conc_model.mean)

            polr_bound_species_type = model.species_types.get_or_create(
                id='{}_bound_{}'.format(polr_complex.id, gene.id),
                name='{} bound {}'.format(polr_complex.name, gene.name),
                type=onto['WC:pseudo_species'],                
                )
            polr_bound_species_type.structure = wc_lang.ChemicalStructure(
                    empirical_formula = polr_complex.structure.empirical_formula,
                    molecular_weight = polr_complex.structure.molecular_weight,
                    charge = polr_complex.structure.charge)
            polr_bound_species = model.species.get_or_create(
                species_type=polr_bound_species_type, compartment=transcription_compartment)
            polr_bound_species.id = polr_bound_species.gen_id()
            self._elongation_modifier[rna_kb.id] = polr_bound_species
            self._gene_bound_polr[rna_pol_pair[rna_kb.id]].append(polr_bound_species)

            conc_model = model.distribution_init_concentrations.create(
                species=polr_bound_species,
                units=unit_registry.parse_units('molecule'),
                )
            conc_model.id = conc_model.gen_id()

            init_reaction = model.reactions.create(
                submodel=self.submodel, id='transcription_initiation_' + rna_kb.id,
                name='transcription initiation of ' + rna_kb.name,
                reversible=False, comments='Set to irreversible to model only the net flux')
            
            init_reaction.participants.append(
                polr_bound_non_specific_species.species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(
                polr_binding_site_species.species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(
                polr_bound_species.species_coefficients.get_or_create(
                coefficient=1))
            init_reaction.participants.append(
                polr_non_specific_binding_site_species.species_coefficients.get_or_create(
                coefficient=1))

            # Add ATP hydrolysis requirement for DNA melting and promoter escape by RNA polymerase II
            if 'RNA Polymerase II' in rna_pol_pair[rna_kb.id]:
                init_reaction.participants.append(metabolites['atp'][
                    transcription_compartment.id].species_coefficients.get_or_create(
                    coefficient=-2))
                init_reaction.participants.append(metabolites['h2o'][
                    transcription_compartment.id].species_coefficients.get_or_create(
                    coefficient=-2))
                init_reaction.participants.append(metabolites['adp'][
                    transcription_compartment.id].species_coefficients.get_or_create(
                    coefficient=2))
                init_reaction.participants.append(metabolites['pi'][
                    transcription_compartment.id].species_coefficients.get_or_create(
                    coefficient=2))
                init_reaction.participants.append(metabolites['h'][
                    transcription_compartment.id].species_coefficients.get_or_create(
                    coefficient=2))

            # Create elongation reaction
            rna_model = model.species_types.get_one(id=rna_kb.id).species[0]
            reaction = model.reactions.get_or_create(
                submodel=self.submodel, id='transcription_elongation_' + rna_kb.id,
                name='transcription elongation of ' + rna_kb.name,
                reversible=False, comments='Lumped reaction')            

            if rna_kb.gene.strand == wc_kb.core.PolymerStrand.positive:
                pre_rna_seq = gene_seq.transcribe()
            else:
                pre_rna_seq = gene_seq.reverse_complement().transcribe()  
            
            if rna_kb.id in gvar.transcript_ntp_usage:
                ntp_count = gvar.transcript_ntp_usage[rna_kb.id]
            else:
                seq = rna_kb.get_seq()
                ntp_count = gvar.transcript_ntp_usage[rna_kb.id] = {
                    'A': seq.upper().count('A'),
                    'C': seq.upper().count('C'),
                    'G': seq.upper().count('G'),
                    'U': seq.upper().count('U'),
                    'len': len(seq)
                    }

            # Adding participants to LHS
            reaction.participants.append(
                polr_bound_species.species_coefficients.get_or_create(
                coefficient=-1))
            reaction.participants.append(metabolites['atp'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=-pre_rna_seq.upper().count('A')))
            reaction.participants.append(metabolites['ctp'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=-pre_rna_seq.upper().count('C')))
            reaction.participants.append(metabolites['gtp'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=-pre_rna_seq.upper().count('G')))
            reaction.participants.append(metabolites['utp'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=-pre_rna_seq.upper().count('U')))
            reaction.participants.append(metabolites['h2o'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=-(len(pre_rna_seq)-pre_rna_seq.upper().count('N')
                    -ntp_count['len']+1)))
            
            # Adding participants to RHS
            reaction.participants.append(
                rna_model.species_coefficients.get_or_create(
                coefficient=1))
            reaction.participants.append(metabolites['ppi'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=len(pre_rna_seq)-pre_rna_seq.upper().count('N')))
            reaction.participants.append(metabolites['amp'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=pre_rna_seq.upper().count('A')-ntp_count['A']))
            reaction.participants.append(metabolites['cmp'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=pre_rna_seq.upper().count('C')-ntp_count['C']))
            reaction.participants.append(metabolites['gmp'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=pre_rna_seq.upper().count('G')-ntp_count['G']))
            reaction.participants.append(metabolites['ump'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=pre_rna_seq.upper().count('U')-ntp_count['U']))
            reaction.participants.append(metabolites['h'][
                transcription_compartment.id].species_coefficients.get_or_create(
                coefficient=len(pre_rna_seq)-pre_rna_seq.upper().count('N')
                    -ntp_count['len']+1))
            reaction.participants.append(
                polr_complex_species.species_coefficients.get_or_create(
                coefficient=1))
            reaction.participants.append(
                polr_binding_site_species.species_coefficients.get_or_create(
                coefficient=1))

            if rna_kb.type==wc_kb.eukaryote.TranscriptType.mRna:
                ribo_binding_site_st = model.species_types.get_or_create(
                    id='{}_ribosome_binding_site'.format(rna_kb.id),
                    name='ribosome binding site of {}'.format(rna_kb.name),
                    type=onto['WC:pseudo_species'],
                    )
                ribo_binding_site_st.structure = wc_lang.ChemicalStructure(
                    empirical_formula = EmpiricalFormula(),
                    molecular_weight = 0.,
                    charge = 0)
                ribo_binding_site_species = model.species.get_or_create(
                    species_type=ribo_binding_site_st, compartment=translation_compartment)
                ribo_binding_site_species.id = ribo_binding_site_species.gen_id()

                site_per_rna = math.floor(ntp_count['len']/ribosome_occupancy_width) + 1
                reaction.participants.append(
                    ribo_binding_site_species.species_coefficients.get_or_create(
                    coefficient=site_per_rna))

                rna_init_conc = model.distribution_init_concentrations.get_one(
                    species=rna_model).mean
                conc_model = model.distribution_init_concentrations.create(
                    species=ribo_binding_site_species,
                    mean=site_per_rna * rna_init_conc,
                    units=unit_registry.parse_units('molecule'),
                    comments='Set to mRNA length divided by {} bp to allow '
                        'queueing of ribosome during translation'.format(ribosome_occupancy_width),
                    references=[ref_ribo_width]    
                    )
                conc_model.id = conc_model.gen_id()

            init_el_rxn_no += 1

        print('{} reactions each for initiation and elongation have been generated'.format(
            init_el_rxn_no))    

    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """

        model = self.model        
        cell = self.knowledge_base.cell
        nucleus = model.compartments.get_one(id='n')  
        mitochondrion = model.compartments.get_one(id='m')

        ref_polr_width = model.references.get_one(
            title='Structure and mechanism of the RNA Polymerase II transcription machinery')

        ref_model = wc_lang.Reference(
            model=model,
            title='Transcriptional regulation by the numbers: models',
            author='Lacramioara Bintu, Nicolas E. Buchler, Hernan G. Garcia, '
                'Ulrich Gerland, Terence Hwa, Jane Kondev, Rob Phillips',        
            year=2005,
            type=onto['WC:article'],
            publication='Current Opinion in Genetics and Development',        
            volume='15',
            pages='116-124'
            )
        ref_model.id = 'ref_'+str(len(model.references))

        ref_kd = wc_lang.Reference(
            model=model,
            title='Macromolecular crowding as a regulator of gene transcription',
            author='Hiroaki Matsuda, Gregory Garbes Putzel, Vadim Backman, Igal Szleifer',        
            year=2014,
            type=onto['WC:article'],
            publication='Biophysical Journal',        
            volume='106',
            pages='1801-1810'
            )
        ref_kd.id = 'ref_'+str(len(model.references))

        Kd_non_specific_polr = model.parameters.create(
            id='K_d_non_specific_polr',
            type=None,
            value=1e-03,
            units=unit_registry.parse_units('M'),
            references=[ref_kd],
            comments='Value taken from the estimation used in the reference'
            )
        
        Kd_specific_polr = model.parameters.create(
            id='K_d_specific_polr',
            type=None,
            value=1e-09,
            units=unit_registry.parse_units('M'),
            references=[ref_kd],
            comments='Value taken from the estimation used in the reference'
            )

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
            comments='Boolean switch for determining if binding site is still available'
            )                         
                
        # Generate rate law for binding of RNA polymerase to non-specific site       
        rna_pol_pair = self.options.get('rna_pol_pair')
        self._polr_pool = {}
        for polr in set(rna_pol_pair.values()):            
            transcription_compartment = mitochondrion if 'mito' in polr else nucleus
            ns_binding_reaction = model.reactions.get_one(
                name='non-specific binding of {} in {}'.format(polr, transcription_compartment.name))

            polr_complex = model.species_types.get_one(name=polr)
            polr_complex_species = model.species.get_one(
                species_type=polr_complex, compartment=transcription_compartment)

            non_specific_binding_constant = model.parameters.create(
                id='k_non_specific_binding_{}'.format(polr_complex.id),
                type=None,
                units=unit_registry.parse_units('molecule^-1 s^-1')
                )

            expression, error = wc_lang.RateLawExpression.deserialize(
                '{} * {}'.format(non_specific_binding_constant.id, polr_complex_species.id), {
                wc_lang.Species: {polr_complex_species.id: polr_complex_species},
                wc_lang.Parameter: {non_specific_binding_constant.id: non_specific_binding_constant},            
                })
            assert error is None, str(error)

            ns_binding_rate_law = model.rate_laws.create(
                direction=wc_lang.RateLawDirection.forward,
                type=None,
                expression=expression,
                reaction=ns_binding_reaction,
                )
            ns_binding_rate_law.id = ns_binding_rate_law.gen_id()

            # Create observable for total RNA polymerase
            polr_bound_non_specific_species = model.species.get_one(
                species_type=model.species_types.get_one(
                    id='{}_bound_non_specific_site'.format(polr_complex.id)), 
                compartment=transcription_compartment)
            
            self._polr_pool[polr] = {i.id: i for i in self._gene_bound_polr[polr]}
            self._polr_pool[polr][polr_complex_species.id] = polr_complex_species        
            self._polr_pool[polr][polr_bound_non_specific_species.id] = polr_bound_non_specific_species 
            
            chunked_data = [[k, v] for k, v in self._polr_pool[polr].items()]
            size = 800
            chunked_data = [chunked_data[i * size:(i + 1) * size] for i in range((len(chunked_data) + size - 1) // size )]
            
            all_subtotal_obs = {}
            for ind, chunk in enumerate(chunked_data):
                chunked_dict = {k:v for k,v in chunk}              
                polr_subtotal_exp, error = wc_lang.ObservableExpression.deserialize(
                    ' + '.join(chunked_dict.keys()),
                    {wc_lang.Species: chunked_dict})            
                assert error is None, str(error)                
                polr_subtotal_obs = model.observables.create(
                    id='subtotal_{}_{}_{}'.format(polr_complex_species.species_type.id, transcription_compartment.id, ind+1), 
                    name='subtotal {} of {} in {}'.format(ind+1, polr, transcription_compartment.name), 
                    units=unit_registry.parse_units('molecule'), 
                    expression=polr_subtotal_exp)
                all_subtotal_obs[polr_subtotal_obs.id] = polr_subtotal_obs
                    
            polr_obs_exp, error = wc_lang.ObservableExpression.deserialize(
                ' + '.join(all_subtotal_obs.keys()),
                {wc_lang.Observable: all_subtotal_obs})            
            assert error is None, str(error)
            polr_obs = model.observables.create(
                id='total_{}_{}'.format(polr_complex_species.species_type.id, transcription_compartment.id), 
                name='total {} in {}'.format(polr, transcription_compartment.name), 
                units=unit_registry.parse_units('molecule'), 
                expression=polr_obs_exp)

            specific_binding_constant = model.parameters.create(
                id='k_specific_binding_{}'.format(polr_complex.id),
                type=None,
                units=unit_registry.parse_units('molecule^-2 s^-1')
                )                    
        
        # Generate response function for initiation, elongation and termination factors for each RNA polymerase
        init_factors = self.options.get('init_factors')
        elongation_termination_factors = self.options.get('elongation_termination_factors')
        elongation_negative_factors = self.options.get('elongation_negative_factors')
        beta = self.options.get('beta')
        
        # Generate response function for each transcription initiation factor group
        init_factor_functions = {}
        for rnap, factors in init_factors.items():
            init_factor_functions[rnap] = {}
            compartment = mitochondrion if 'polm' in rnap else nucleus
            n = 1
            for factor in factors:
                factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                    model, beta, 'transcription_init_{}'.format(rnap[:4]), 'transcription_init_{}'.format(rnap[:4]), 
                    compartment, [factor])
                
                objects = {
                            wc_lang.Species: all_species,
                            wc_lang.Parameter: all_parameters,
                            wc_lang.Observable: all_observables,
                            wc_lang.Function: all_volumes,            
                            }                                
                        
                expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                assert error is None, str(error)

                init_factor_functions[rnap][','.join(factor)] = {
                    'function': model.functions.create(     
                        id='transcription_init_factor_function_{}_{}'.format(rnap, n),               
                        name='response function for transcription initiation factor {} for {}'.format(n, rnap),
                        expression=expression,
                        units=unit_registry.parse_units(''),
                        ),
                    'objects': objects}
                n += 1
        
        # Generate response function for each transcription elongation factor group
        elongation_termination_factor_functions = {}
        for rnap, factors in elongation_termination_factors.items():
            elongation_termination_factor_functions[rnap] = {}
            compartment = mitochondrion if 'polm' in rnap else nucleus
            n = 1
            for factor in factors:
                factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                    model, beta, 'transcription_el_{}'.format(rnap[:4]), 'transcription_el_{}'.format(rnap[:4]), 
                    compartment, [factor])
                
                objects = {
                            wc_lang.Species: all_species,
                            wc_lang.Parameter: all_parameters,
                            wc_lang.Observable: all_observables,
                            wc_lang.Function: all_volumes,            
                            }                                
                        
                expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                assert error is None, str(error)

                elongation_termination_factor_functions[rnap][','.join(factor)] = {
                    'function': model.functions.create(     
                        id='transcription_el_factor_function_{}_{}'.format(rnap, n),               
                        name='response function for transcription elongation factor {} for {}'.format(n, rnap),
                        expression=expression,
                        units=unit_registry.parse_units(''),
                        ),
                    'objects': objects}
                n += 1

        # Generate response function for each transcription elongation negative factor group (only RNA Pol II has this factor)
        elongation_negative_factor_functions = {}
        for rnap, factors in elongation_negative_factors.items():
            elongation_negative_factor_functions[rnap] = {}
            compartment = mitochondrion if 'polm' in rnap else nucleus
            n = 1
            for factor in factors:
                factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                    model, beta, 'transcription_neg_{}'.format(rnap[:4]), 'transcription_neg_{}'.format(rnap[:4]), 
                    compartment, [factor])
                
                objects = {
                            wc_lang.Species: all_species,
                            wc_lang.Parameter: all_parameters,
                            wc_lang.Observable: all_observables,
                            wc_lang.Function: all_volumes,            
                            }                                
                        
                expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                assert error is None, str(error)

                elongation_negative_factor_functions[rnap][','.join(factor)] = {
                    'function': model.functions.create(     
                        id='transcription_neg_factor_function_{}_{}'.format(rnap, n),               
                        name='response function for transcription elongation negative factor {} for {}'.format(n, rnap),
                        expression=expression,
                        units=unit_registry.parse_units(''),
                        ),
                    'objects': objects}
                n += 1

        # Generate rate laws for initiation and elongation & termination
        rna_init_factors = self.options.get('rna_init_factors')
        rna_elongation_termination_factors = self.options.get('rna_elongation_termination_factors')
        rna_elongation_negative_factors = self.options.get('rna_elongation_negative_factors')
        polr_occupancy_width = self.options.get('polr_occupancy_width')
        p_function_exprs = {}
        self._gene_p_function_map = {}
        rate_law_no = 0                    
        rnas_kb = cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType)
        for rna_kb in rnas_kb:

            rna_kb_compartment_id = rna_kb.species[0].compartment.id
            if rna_kb_compartment_id == 'c':
                no_of_binding_sites = model.parameters.get_or_create(
                    id='total_nuclear_genome_binding',
                    type=None,
                    value=self._nuclear_max_binding_sites,
                    units=unit_registry.parse_units('molecule'),
                    references=[ref_polr_width],
                    comments='Set to genome length divided by {} bp'.format(polr_occupancy_width)
                    )
                transcription_compartment = nucleus      
            else:
                no_of_binding_sites = model.parameters.get_or_create(
                    id='total_mitochondrial_genome_binding',
                    type=None,
                    value=self._mitochondrial_max_binding_sites,
                    units=unit_registry.parse_units('molecule'),
                    references=[ref_polr_width],
                    comments='Set to genome length divided by {} bp'.format(polr_occupancy_width)
                    )
                transcription_compartment = mitochondrion

            # Assign transcriptional regulation
            reg_species = {}            
            reg_parameters = {}
            reg_functions = {}

            reg_parameters[Kd_specific_polr.id] = Kd_specific_polr
            reg_parameters[Kd_non_specific_polr.id] = Kd_non_specific_polr    
            reg_parameters[no_of_binding_sites.id] = no_of_binding_sites        
            
            F_regs = []
            reaction_id = 'transcription_initiation_' + rna_kb.id
            for reg in rna_kb.gene.regulatory_modules:                                
                for tf in reg.transcription_factor_regulation: 
                    
                    tf_model = model.species.get_one(
                        species_type=model.species_types.get_one(id=tf.transcription_factor.id), 
                        compartment=transcription_compartment)                   
                    
                    if tf_model and tf.direction == wc_kb.eukaryote.RegulatoryDirection.activation:                        
                        F_act, species_act, param_act, func_act = utils.simple_activator(
                            model, reaction_id, tf_model)
                        F_regs.append(F_act)
                        reg_species.update(species_act)
                        reg_parameters.update(param_act)
                        reg_functions.update(func_act)
                        
                    elif tf_model and tf.direction == wc_kb.eukaryote.RegulatoryDirection.repression:
                        F_rep, species_rep, param_rep, func_rep = utils.simple_repressor(
                            model, reaction_id, tf_model) 
                        F_regs.append(F_rep)
                        reg_species.update(species_rep)
                        reg_parameters.update(param_rep)
                        reg_functions.update(func_rep)
            
            F_reg_N = ' * '.join(F_regs)

            # Generate rate law for initiation
            polr_bound_non_specific_species = model.species.get_one(
                species_type=model.species_types.get_one(
                    id='{}_bound_non_specific_site'.format(
                        self._initiation_polr_species[rna_kb.id].species_type.id)), 
                compartment=transcription_compartment)
            reg_species[polr_bound_non_specific_species.id] = polr_bound_non_specific_species

            polr_complex_species = model.species.get_one(
                species_type=model.species_types.get_one(name=rna_pol_pair[rna_kb.id]), 
                compartment=transcription_compartment)
            polr_obs = model.observables.get_one(
                id='total_{}_{}'.format(polr_complex_species.species_type.id, transcription_compartment.id))

            p_bound = '1 / (1 + {} / ({} * {}) * exp(log({} / {})))'.format(
                no_of_binding_sites.id, 
                polr_obs.id, 
                F_reg_N if F_reg_N else 1,
                Kd_specific_polr.id,
                Kd_non_specific_polr.id
                )
            if p_bound in p_function_exprs:
                p_bound_function = p_function_exprs[p_bound]
                self._gene_p_function_map[rna_kb.id] = p_bound_function
            else:    
                p_bound_expression, error = wc_lang.FunctionExpression.deserialize(p_bound, {
                    wc_lang.Species: reg_species,
                    wc_lang.Parameter: reg_parameters,            
                    wc_lang.Function: reg_functions,
                    wc_lang.Observable: {polr_obs.id: polr_obs}
                    })
                assert error is None, str(error)
                
                p_bound_function = model.functions.create(                    
                    name='probability of RNAP binding {}'.format(len(p_function_exprs)),
                    expression=p_bound_expression,
                    references=[ref_model],
                    units=unit_registry.parse_units(''),
                    )
                p_bound_function.id = 'p_bound_{}'.format(len(p_function_exprs)+1)
                p_function_exprs[p_bound] = p_bound_function
                self._gene_p_function_map[rna_kb.id] = p_bound_function
            
            specific_binding_constant = model.parameters.get_one(
                id='k_specific_binding_{}'.format(polr_complex_species.species_type.id))
            reg_parameters[specific_binding_constant.id] = specific_binding_constant

            objects = {
                wc_lang.Species: {},
                wc_lang.Parameter: {},
                wc_lang.Observable: {},
                wc_lang.Function: {},            
                }
            expression_terms = []
            for factor in init_factors[rna_init_factors[rna_kb.id]]:
                factor_details = init_factor_functions[rna_init_factors[rna_kb.id]][','.join(factor)]
                expression_terms.append(factor_details['function'].id)
                for cl, dictionary in objects.items():
                    dictionary.update(factor_details['objects'][cl])
                objects[wc_lang.Function][factor_details['function'].id] = factor_details['function']

            expression = '{} * {} * {} * max(min({} , {}) , {}) * {} * 2**{}'.format(
                p_bound_function.id,
                specific_binding_constant.id,
                polr_bound_non_specific_species.id,
                self._allowable_queue_len[rna_kb.id][0].id,
                max_bool.id,
                min_bool.id,
                ' * '.join(expression_terms),
                len(expression_terms),                
                )
            reg_species[self._allowable_queue_len[rna_kb.id][0].id] = self._allowable_queue_len[rna_kb.id][0]
            reg_parameters[max_bool.id] = max_bool
            reg_parameters[min_bool.id] = min_bool

            objects[wc_lang.Species].update(reg_species)
            objects[wc_lang.Parameter].update(reg_parameters)
            objects[wc_lang.Function].update({p_bound_function.id: p_bound_function})
            
            init_rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, objects)
            assert error is None, str(error)

            init_rate_law = model.rate_laws.create(
                direction=wc_lang.RateLawDirection.forward,
                type=None,
                expression=init_rate_law_expression,
                reaction=model.reactions.get_one(id='transcription_initiation_' + rna_kb.id),
                units=unit_registry.parse_units('s^-1'),                
                )
            init_rate_law.id = init_rate_law.gen_id()                   

            # Generate rate law for the lumped reaction of elongation & termination
            elongation_reaction = model.reactions.get_one(id='transcription_elongation_' + rna_kb.id)

            objects = {
                wc_lang.Species: {},
                wc_lang.Parameter: {},
                wc_lang.Observable: {},
                wc_lang.Function: {},            
                }
            expression_terms = []
            for factor in elongation_termination_factors[rna_elongation_termination_factors[rna_kb.id]]:
                factor_details = elongation_termination_factor_functions[
                    rna_elongation_termination_factors[rna_kb.id]][','.join(factor)]
                expression_terms.append(factor_details['function'].id)
                for cl, dictionary in objects.items():
                    dictionary.update(factor_details['objects'][cl])
                objects[wc_lang.Function][factor_details['function'].id] = factor_details['function']

            if elongation_negative_factors.get(rna_elongation_negative_factors[rna_kb.id]):
                for factor in elongation_negative_factors[rna_elongation_negative_factors[rna_kb.id]]:
                    factor_details = elongation_negative_factor_functions[
                        rna_elongation_negative_factors[rna_kb.id]][','.join(factor)]
                    expression_terms.append(factor_details['function'].id)
                    for cl, dictionary in objects.items():
                        dictionary.update(factor_details['objects'][cl])
                    objects[wc_lang.Function][factor_details['function'].id] = factor_details['function']    
            
            polr_gene_bound_species = self._elongation_modifier[rna_kb.id]
            objects[wc_lang.Species][polr_gene_bound_species.id] = polr_gene_bound_species

            substrates = [[i.species_type.id] for i in elongation_reaction.get_reactants() 
                if (i.species_type.id!='h2o' and i!=polr_gene_bound_species)]
            expressions, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                model, beta, elongation_reaction.id, 'transcription_elongation', transcription_compartment, substrates)
            expression_terms += expressions
            objects[wc_lang.Species].update(all_species)
            objects[wc_lang.Parameter].update(all_parameters)
            objects[wc_lang.Function].update(all_volumes)
            objects[wc_lang.Observable].update(all_observables)            
            
            k_cat_elongation = model.parameters.create(
                id='k_cat_{}'.format(elongation_reaction.id),
                type=onto['WC:k_cat'],
                units=unit_registry.parse_units('molecule^-1 s^-1'),
                )
            objects[wc_lang.Parameter][k_cat_elongation.id] = k_cat_elongation

            expression = '{} * {} * {} * 2**{}'.format(
                k_cat_elongation.id,
                polr_gene_bound_species.id,
                ' * '.join(expression_terms),
                len(expression_terms),
                )

            el_rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, objects)
            assert error is None, str(error)
            
            rate_law = model.rate_laws.create(
                direction=wc_lang.RateLawDirection.forward,
                type=None,
                expression=el_rate_law_expression,
                reaction=elongation_reaction,
                )
            rate_law.id = rate_law.gen_id()

            rate_law_no += 1

        print('{} rate laws for initiation and elongation have been generated'.format(rate_law_no))              
        
    def calibrate_submodel(self):
        """ Calibrate the submodel using data in the KB """
        
        model = self.model        
        cell = self.knowledge_base.cell
        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')

        beta = self.options.get('beta')
        beta_activator = self.options.get('beta_activator')
        beta_repressor = self.options.get('beta_repressor')
        activator_effect = self.options.get('activator_effect')
        rna_pol_pair = self.options.get('rna_pol_pair')

        Avogadro = model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))

        mean_doubling_time = model.parameters.get_one(id='mean_doubling_time').value

        average_rate = {}
        p_bound = {}
        rnas_kb = cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType)
        for rna_kb in rnas_kb:            
        
            transcription_compartment = nucleus if rna_kb.species[0].compartment.id == 'c' else mitochondrion 
            
            # Estimate the average rate of transcription
            rna_product = model.species_types.get_one(id=rna_kb.id).species[0]
            
            half_life = rna_kb.properties.get_one(property='half-life').get_value()
            mean_concentration = rna_product.distribution_init_concentration.mean         

            average_rate[rna_kb.id] = utils.calc_avg_syn_rate(
                mean_concentration, half_life, mean_doubling_time)

            # Estimate the average probability of RNA polymerase binding            
            init_reg_species_count = {}                
            init_reaction = model.reactions.get_one(id='transcription_initiation_' + rna_kb.id)
            for param in init_reaction.rate_laws[0].expression.functions[0].expression.parameters:
                
                if 'Kr_' in param.id:
                    repressor_species = model.species.get_one(
                        id='{}[{}]'.format(param.id.split('_')[-1], transcription_compartment.id))
                    init_reg_species_count[repressor_species.id] = \
                        repressor_species.distribution_init_concentration.mean
                    if repressor_species.distribution_init_concentration.mean:    
                        param.value = beta_repressor * repressor_species.distribution_init_concentration.mean \
                            / Avogadro.value / repressor_species.compartment.init_volume.mean
                        param.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                            beta_repressor, repressor_species.species_type.id, repressor_species.compartment.name)
                    else:
                        param.value = 1e-05
                        param.comments = 'The value was assigned to 1e-05 because the concentration of {} in {} was zero'.format(
                            repressor_species.species_type.id, repressor_species.compartment.name)
                        
                elif 'Ka_' in param.id:
                    activator_species = model.species.get_one(
                        id='{}[{}]'.format(param.id.split('_')[-1], transcription_compartment.id))
                    init_reg_species_count[activator_species.id] = \
                        activator_species.distribution_init_concentration.mean
                    if activator_species.distribution_init_concentration.mean:    
                        param.value = beta_activator * activator_species.distribution_init_concentration.mean \
                            / Avogadro.value / activator_species.compartment.init_volume.mean
                        param.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                            beta_activator, activator_species.species_type.id, activator_species.compartment.name) 
                    else:
                        param.value = 1e-05
                        param.comments = 'The value was assigned to 1e-05 because the concentration of {} in {} was zero'.format(
                            activator_species.species_type.id, activator_species.compartment.name)                      
                
                elif 'f_' in param.id:
                    param.value = activator_effect

            total_polr = self._total_polr[rna_pol_pair[rna_kb.id]]
            no_of_polr_pool = len(self._polr_pool[rna_pol_pair[rna_kb.id]])
            for i in self._polr_pool[rna_pol_pair[rna_kb.id]]:
                init_reg_species_count[i] = total_polr / no_of_polr_pool            
            
            p_bound_function = self._gene_p_function_map[rna_kb.id]
            p_bound_value = p_bound_function.expression._parsed_expression.eval({
                wc_lang.Species: init_reg_species_count,
                wc_lang.Compartment: {
                    transcription_compartment.id: transcription_compartment.init_volume.mean * \
                        transcription_compartment.init_density.value},
                })
            p_bound[rna_kb.id] = p_bound_value
        
        # Calibrate binding constants
        polr_rna_pair = collections.defaultdict(list)
        for rna_id, polr in rna_pol_pair.items():
            polr_rna_pair[polr].append(rna_id)
        
        total_p_bound = {}
        total_gene_bound = {}
        for polr, rnas in polr_rna_pair.items():
            
            transcription_compartment = mitochondrion if 'mito' in polr else nucleus

            polr_complex = model.species_types.get_one(name=polr)
            polr_complex_species = model.species.get_one(
                species_type=polr_complex, compartment=transcription_compartment)
            polr_free_conc = model.distribution_init_concentrations.get_one(
                species=polr_complex_species).mean

            polr_ns_bound = model.species_types.get_one(
                    id='{}_bound_non_specific_site'.format(polr_complex.id))
            polr_ns_bound_species = model.species.get_one(
                species_type=polr_ns_bound, compartment=transcription_compartment)
            polr_ns_bound_conc = model.distribution_init_concentrations.get_one(
                species=polr_ns_bound_species).mean

            total_gene_bound[polr] = self._total_polr[polr] - polr_free_conc - polr_ns_bound_conc

            total_polr_usage_rate = 0
            total_p_bound[polr] = 0
            for rna_id in rnas:
                total_polr_usage_rate += average_rate[rna_id]
                total_p_bound[polr] += p_bound[rna_id]           
            
            non_specific_binding_constant = model.parameters.get_one(
                id='k_non_specific_binding_{}'.format(polr_complex.id))
            non_specific_binding_constant.value = total_polr_usage_rate / polr_free_conc

            specific_binding_constant = model.parameters.get_one(
                id='k_specific_binding_{}'.format(polr_complex.id))
            specific_binding_constant.value = total_polr_usage_rate / \
                (polr_ns_bound_conc * total_p_bound[polr])
            
        # Calibrate the reaction constant of lumped elongation and termination                         
        undetermined_model_kcat = []
        determined_kcat = []
        for rna_kb in rnas_kb: 

            transcription_compartment = nucleus if rna_kb.species[0].compartment.id == 'c' else mitochondrion

            polr_gene_bound_conc = min(self._allowable_queue_len[rna_kb.id][1], 
                p_bound[rna_kb.id] / total_p_bound[rna_pol_pair[rna_kb.id]] * \
                total_gene_bound[rna_pol_pair[rna_kb.id]])
            
            polr_gene_bound_species = self._elongation_modifier[rna_kb.id]            
            model.distribution_init_concentrations.get_one(
                species=polr_gene_bound_species).mean = polr_gene_bound_conc
            
            init_species_counts = {}
            reaction = model.reactions.get_one(id='transcription_elongation_' + rna_kb.id)
            for species in reaction.rate_laws[0].expression.species:
                init_species_counts[species.id] = species.distribution_init_concentration.mean
                if model.parameters.get(id='K_m_{}_{}'.format(reaction.id, species.species_type.id)):
                    model_Km = model.parameters.get_one(
                        id='K_m_{}_{}'.format(reaction.id, species.species_type.id))
                    model_Km.value = beta * species.distribution_init_concentration.mean \
                        / Avogadro.value / species.compartment.init_volume.mean
                    model_Km.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                        beta, species.species_type.id, species.compartment.name)

            for func in reaction.rate_laws[0].expression.functions:
                for species in func.expression.species:
                    init_species_counts[species.id] = species.distribution_init_concentration.mean
                for obs in func.expression.observables:
                    for species in obs.expression.species:    
                        init_species_counts[species.id] = species.distribution_init_concentration.mean            
            
            model_kcat = model.parameters.get_one(id='k_cat_{}'.format(reaction.id))

            if polr_gene_bound_conc and average_rate[rna_kb.id]:
                model_kcat.value = 1.
                model_kcat.value = average_rate[rna_kb.id] / \
                    reaction.rate_laws[0].expression._parsed_expression.eval({
                        wc_lang.Species: init_species_counts,
                        wc_lang.Compartment: {
                            transcription_compartment.id: transcription_compartment.init_volume.mean * \
                                transcription_compartment.init_density.value},
                    })
                determined_kcat.append(model_kcat.value)
            else:
                undetermined_model_kcat.append(model_kcat)
        
        median_kcat = numpy.median(determined_kcat)
        for model_kcat in undetermined_model_kcat:
            model_kcat.value = median_kcat
            model_kcat.comments = 'Set to the median value because it could not be determined from data'

        print('Transcription submodel has been generated')        
