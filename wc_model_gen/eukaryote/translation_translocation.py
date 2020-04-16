""" Generator for translation, protein folding and translocation submodel for eukaryotes
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-06-14
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import wc_model_gen.global_vars as gvar
import wc_model_gen.utils as utils
import Bio.Alphabet
import Bio.Seq
import math
import numpy
import scipy.constants
import wc_kb
import wc_lang
import wc_model_gen

ANTICODON_CODON_RECOGNITION_RULES = {
    'GAA': ['TTT', 'TTC'],
    'TAA': ['TTA', 'TTG'],
    'CAA': ['TTG'],
    'AGA': ['TCT', 'TCC', 'TCA'],
    'GGA': ['TCT', 'TCC'],
    'TGA': ['TCA', 'TCG'],
    'CGA': ['TCG'],
    'GTA': ['TAT', 'TAC'],
    'GCA': ['TGT', 'TGC'],
    'CCA': ['TGG'],
    'AAG': ['CTT', 'CTC', 'CTA'],
    'GAG': ['CTT', 'CTC'],
    'TAG': ['CTA', 'CTG'],
    'CAG': ['CTG'],
    'AGG': ['CCT', 'CCC', 'CCA'],
    'GGG': ['CCT', 'CCC'],
    'TGG': ['CCA', 'CCG'],
    'CGG': ['CCG'],
    'GTG': ['CAT', 'CAC'],
    'TTG': ['CAA', 'CAG'],
    'CTG': ['CAG'],
    'ACG': ['CGT', 'CGC', 'CGA'],
    'GCG': ['CGT', 'CGC'],
    'TCG': ['CGA', 'CGG'],
    'CCG': ['CGG'],
    'AAT': ['ATT', 'ATC', 'ATA'],
    'GAT': ['ATT', 'ATC', 'ATA'],
    'TAT': ['ATA'],
    'CAT': ['ATG'],
    'AGT': ['ACT', 'ACC', 'ACA'],
    'GGT': ['ACT', 'ACC'],
    'TGT': ['ACA', 'ACG'],
    'CGT': ['ACG'],
    'GTT': ['AAT', 'AAC'],
    'TTT': ['AAA', 'AAG'],
    'CTT': ['AAG'],
    'GCT': ['AGT', 'AGC'],
    'TCT': ['AGA', 'AGG'],
    'CCT': ['AGG'],
    'AAC': ['GTT', 'GTC', 'GTA'],
    'GAC': ['GTT', 'GTC'],
    'TAC': ['GTA', 'GTG'],
    'CAC': ['GTG'],
    'AGC': ['GCT', 'GCC', 'GCA'],
    'GGC': ['GCT', 'GCC'],
    'TGC': ['GCA', 'GCG'],
    'CGC': ['GCG'],
    'GTC': ['GAT', 'GAC'],
    'TTC': ['GAA', 'GAG'],
    'CTC': ['GAG'],
    'ACC': ['GGT', 'GGC', 'GGA'],
    'GCC': ['GGT', 'GGC'],
    'TCC': ['GGA', 'GGG'],
    'CCC': ['GGG'],
    'TCA': ['TGA'], #selenocysteine
    'AAA': ['TTT'], #natural pairing but unlikely according to the rule
    'ATA': ['TAT'], #natural pairing but unlikely according to the rule
    'ATG': ['CAT'], #natural pairing but unlikely according to the rule
    'ATT': ['AAT'], #natural pairing but unlikely according to the rule
    'ACT': ['AGT'], #natural pairing but unlikely according to the rule
    'ATC': ['GAT'], #natural pairing but unlikely according to the rule        
}
UNCONDITIONAL_STOP_CODON = ['TAA', 'TAG']
CONDITIONAL_STOP_CODON = ['TGA']


class TranslationTranslocationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for translation, protein folding and translocation submodel

        Translation, protein folding and translocation processes are 
        modeled as three reaction steps in this submodel:

        1. Translation initiation where ribosomes and methionine (or other start amino acid) 
           bind to the mRNA. For nuclear mRNAs, transport from the nucleus to the cytoplasm 
           are lumped with this reaction. The energetic of met-tRNA charging is included;
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
        * mitochondrial_exosome (:obj:`str`): the name of exosome complex that degrades RNAs in 
            the mitochondria                    
        * amino_acid_id_conversion (:obj:`dict`): a dictionary with amino acid standard ids
            as keys and amino acid metabolite ids as values     
        * codon_table (:obj:`dict`, optional): a dictionary with protein id as key and 
            NCBI identifier for translation table as value, the default is 1 (standard table) 
            for all protein
        * cds (:obj:`bool`, optional): True indicates the sequences of protein are complete CDS,
            the default is True        
        * beta (:obj:`float`, optional): ratio of Michaelis-Menten constant to substrate 
            concentration (Km/[S]) for use when estimating Km values, the default value is 1
        * polysome_fraction (:obj:`dict`): a dictionary with mRNA ids as keys and
            fraction of total cellular ribosomes the mRNA is bound to
        * mitochondrial_cytosolic_trna_partition (:obj:`float`, optional): fraction of cellular 
            tRNA that would be imported into the mitochondrial for codons not covered by the 
            mitochondrial tRNAs, the default value is 0.01
        * selenoproteome (:obj:`list`, optional): list of IDs of genes that translate into 
            selenoproteins, default is an empty list                    
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

        if 'mitochondrial_exosome' not in options:
            raise ValueError('The name of mitochondrial exosome has not been provided')
        else:    
            mitochondrial_exosome = options['mitochondrial_exosome']

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

        mitochondrial_cytosolic_trna_partition = options.get('mitochondrial_cytosolic_trna_partition', 0.01)
        assert(0. <= mitochondrial_cytosolic_trna_partition <= 1.) 
        options['mitochondrial_cytosolic_trna_partition'] = mitochondrial_cytosolic_trna_partition

        selenoproteome = options.get('selenoproteome', [])
        options['selenoproteome'] = selenoproteome    

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        
        cytoplasmic_ribosome = self.options.get('cytoplasmic_ribosome')
        mitochondrial_ribosome = self.options.get('mitochondrial_ribosome')
        amino_acid_id_conversion = self.options.get('amino_acid_id_conversion')
        codon_table = self.options['codon_table']
        cds = self.options['cds']
        selenoproteome = self.options['selenoproteome']        
        
        cytosol = model.compartments.get_one(id='c')
        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')
        peroxisome = model.compartments.get_one(id='x')
                           
        # Get metabolite species involved in reaction
        amino_acid_participants = list(amino_acid_id_conversion.values()) 
        other_metabolite_participants = ['atp', 'adp', 'amp', 'gtp', 'gdp', 'pi', 'ppi', 'h2o', 'h', 'selnp']
        metabolites = {}
        for met in amino_acid_participants + other_metabolite_participants:
            met_species_type = model.species_types.get_one(id=met)
            metabolites[met] = {
                'c': met_species_type.species.get_one(compartment=cytosol),
                'm': met_species_type.species.get_one(compartment=mitochondrion)
                }

        self.submodel.framework = wc_ontology['WC:next_reaction_method']

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
            if mrna_kb_compartment_id == 'c':
                translation_compartment = cytosol
                ribosome_complex = model.species_types.get_one(
                    name=cytoplasmic_ribosome).species.get_one(compartment=cytosol)
            else:
                translation_compartment = mitochondrion
                ribosome_complex = model.species_types.get_one(
                    name=mitochondrial_ribosome).species.get_one(compartment=mitochondrion)            

            # Create initiation reaction
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

            aa_content = {}
            if mrna_kb.protein.id in gvar.protein_aa_usage:                                        
                for aa, aa_id in amino_acid_id_conversion.items():
                    if gvar.protein_aa_usage[mrna_kb.protein.id][aa]:
                        aa_content[aa_id] = gvar.protein_aa_usage[mrna_kb.protein.id][aa]
            else:
                gvar.protein_aa_usage[mrna_kb.protein.id] = {i:0 for i in list(amino_acid_id_conversion.keys())}
                if codon_table == 1:
                    codon_id = 1
                else:
                    codon_id = codon_table[mrna_kb.protein.id]
                _, raw_seq, start_codon = mrna_kb.protein.get_seq_and_start_codon(table=codon_id, cds=cds)
                if mrna_kb.gene.id in selenoproteome:
                    processed_seq = raw_seq[:-1] if raw_seq.endswith('*') else raw_seq
                    protein_seq = ''.join(i if i!='*' else 'U' for i in processed_seq)
                else:                                            
                    protein_seq = ''.join(i for i in raw_seq if i!='*')
                for aa in protein_seq:
                    aa_id = amino_acid_id_conversion[aa]
                    if aa_id not in aa_content:
                        aa_content[aa_id] = 1
                        gvar.protein_aa_usage[mrna_kb.protein.id][aa] = 1
                    else:
                        aa_content[aa_id] += 1
                        gvar.protein_aa_usage[mrna_kb.protein.id][aa] += 1
                gvar.protein_aa_usage[mrna_kb.protein.id]['*'] = raw_seq.count('*')
                gvar.protein_aa_usage[mrna_kb.protein.id]['len'] = len(protein_seq)
                gvar.protein_aa_usage[mrna_kb.protein.id]['start_aa'] = protein_seq[0]
                gvar.protein_aa_usage[mrna_kb.protein.id]['start_codon'] = str(start_codon).upper()

            first_aa = model.species_types.get_one(id=amino_acid_id_conversion[
                gvar.protein_aa_usage[mrna_kb.protein.id]['start_aa']])           
            
            ribo_binding_site_species = model.species_types.get_one(
                id='{}_ribosome_binding_site'.format(mrna_kb.id)).species[0]
            self._allowable_queue_len[mrna_kb.id] = (ribo_binding_site_species, 
                ribo_binding_site_species.distribution_init_concentration.mean)

            ribo_bound_species_type = model.species_types.get_or_create(
                id='ribo_bound_{}'.format(mrna_kb.id),
                name='Ribosome bound {}'.format(mrna_kb.name),
                type=wc_ontology['WC:pseudo_species'],                
                )
            ribo_bound_species_type.structure = wc_lang.ChemicalStructure(
                empirical_formula = ribosome_complex.species_type.structure.empirical_formula +\
                    first_aa.structure.empirical_formula,
                molecular_weight = ribosome_complex.species_type.structure.molecular_weight +\
                    first_aa.structure.molecular_weight,
                charge = ribosome_complex.species_type.structure.charge +\
                    first_aa.structure.charge,
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
            # Include 1 ATP hydrolysis for the charging of tRNA-met (or other start amino acid)            
            init_reaction.participants.append(
                ribosome_complex.species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(
                ribo_binding_site_species.species_coefficients.get_or_create(
                coefficient=-1))
            init_reaction.participants.append(first_aa.species.get_one(
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
                model=model, compartment=translation_compartment)
            protein_model.id = protein_model.gen_id()
            
            if not protein_model.distribution_init_concentration:
                conc_model = model.distribution_init_concentrations.create(
                    species=protein_model,
                    mean=0.,
                    units=unit_registry.parse_units('molecule'),
                    comments='Created and set to zero because the protein is translated ' +\
                        'but not localized in this compartment'
                    )
                conc_model.id = conc_model.gen_id()
            
            el_reaction = model.reactions.get_or_create(
                submodel=self.submodel, id='translation_elongation_' + mrna_kb.id,
                name='translation elongation of ' + mrna_kb.name,
                reversible=False, comments='Lumped reaction')                    
            
            aa_content[amino_acid_id_conversion[gvar.protein_aa_usage[mrna_kb.protein.id]['start_aa']]] -= 1

            # Adding participation of amino acids and other additional metabolites for forming selenocysteine
            serine_no = 0
            for aa_met, count in aa_content.items():   
                if count:
                    # To add selenocysteine, seryl-tRNA is formed and phosphorylated before reacting with selenophosphate
                    if aa_met == amino_acid_id_conversion['U']:
                        serine_no = count
                        serine_met = amino_acid_id_conversion['S']
                        serine_species = metabolites[serine_met][translation_compartment.id]
                        serine_coefficient = el_reaction.participants.get_one(
                            species=serine_species)
                        if serine_coefficient:
                            old_coef = serine_coefficient.coefficient
                            el_reaction.participants.remove(serine_coefficient)
                            el_reaction.participants.add(
                                serine_species.species_coefficients.get_or_create(
                                    coefficient=old_coef - count))
                        else:
                            el_reaction.participants.append(
                                serine_species.species_coefficients.get_or_create(
                                    coefficient=-count))
                        el_reaction.participants.append(metabolites['selnp'][
                            translation_compartment.id].species_coefficients.get_or_create(
                                coefficient=-count))                        
                        el_reaction.participants.append(metabolites['adp'][
                            translation_compartment.id].species_coefficients.get_or_create(
                                coefficient=count))
                        el_reaction.participants.append(metabolites['ppi'][
                            translation_compartment.id].species_coefficients.get_or_create(
                                coefficient=count))
                    else:
                        aa_species = metabolites[aa_met][translation_compartment.id]
                        aa_coefficient = el_reaction.participants.get_one(
                            species=aa_species)
                        if aa_coefficient:
                            old_coef = aa_coefficient.coefficient
                            el_reaction.participants.remove(aa_coefficient)
                            el_reaction.participants.add(
                                aa_species.species_coefficients.get_or_create(
                                    coefficient=old_coef - count))
                        else:             
                            el_reaction.participants.append(
                                aa_species.species_coefficients.get_or_create(
                                    coefficient=-count))

            # Adding general participants to LHS
            # Include 1 ATP hydrolysis for each tRNA-aa charging
            # Include 1 GTP hydrolysis for each peptide bond formation
            # Include 1 GTP hydrolysis at termination 
            el_reaction.participants.append(ribo_bound_species.species_coefficients.get_or_create(
                coefficient=-1))
            el_reaction.participants.append(metabolites['gtp'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=-gvar.protein_aa_usage[mrna_kb.protein.id]['len']))
            el_reaction.participants.append(metabolites['atp'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=-(gvar.protein_aa_usage[mrna_kb.protein.id]['len']- 1 + serine_no)))            
            el_reaction.participants.append(metabolites['h2o'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=-((gvar.protein_aa_usage[mrna_kb.protein.id]['len']-1)*2 + 1)))           
            
            # Adding general participants to RHS
            el_reaction.participants.append(ribosome_complex.species_coefficients.get_or_create(
                coefficient=1))
            el_reaction.participants.append(ribo_binding_site_species.species_coefficients.get_or_create(
                coefficient=1))
            el_reaction.participants.append(protein_model.species_coefficients.get_or_create(
                coefficient=1))
            el_reaction.participants.append(metabolites['amp'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=gvar.protein_aa_usage[mrna_kb.protein.id]['len']-1))
            el_reaction.participants.append(metabolites['gdp'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=gvar.protein_aa_usage[mrna_kb.protein.id]['len']))
            el_reaction.participants.append(metabolites['pi'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=(gvar.protein_aa_usage[mrna_kb.protein.id]['len']-1)*3 + 1))            
            el_reaction.participants.append(metabolites['h'][
                translation_compartment.id].species_coefficients.get_or_create(
                    coefficient=gvar.protein_aa_usage[mrna_kb.protein.id]['len'] - serine_no))

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

                if compartment.id=='n':
                    energy_compartment = nucleus # GTP-dependent translocation
                    energy_reactant = 'gtp'
                    energy_product = 'gdp'
                elif compartment.id=='m':
                    energy_compartment = mitochondrion # ATP-dependent translocation
                    energy_reactant = 'atp'
                    energy_product = 'adp'
                elif compartment.id=='x':
                    energy_compartment = peroxisome # ATP-dependent translocation
                    energy_reactant = 'atp'
                    energy_product = 'adp'
                else:
                    energy_compartment = cytosol # GTP-dependent translocation to other organelles and membranes through er
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
        model = self.model
        cell = self.knowledge_base.cell

        amino_acid_id_conversion = self.options.get('amino_acid_id_conversion')
        beta = self.options.get('beta')
        codon_table = self.options['codon_table']
        cds = self.options['cds']
        selenoproteome = self.options['selenoproteome']

        cytoplasmic_ribosome = self.options.get('cytoplasmic_ribosome')
        mitochondrial_ribosome = self.options.get('mitochondrial_ribosome')
        cytoplasmic_initiation_factors = self.options.get('cytoplasmic_initiation_factors')
        mitochondrial_initiation_factors = self.options.get('mitochondrial_initiation_factors')
        cytoplasmic_elongation_factors = self.options.get('cytoplasmic_elongation_factors')
        mitochondrial_elongation_factors = self.options.get('mitochondrial_elongation_factors')
        cytoplasmic_chaperones = self.options.get('cytoplasmic_chaperones')
        mitochondrial_chaperones = self.options.get('mitochondrial_chaperones')
        er_chaperones = self.options.get('er_chaperones')
        
        cytosol = model.compartments.get_one(id='c')
        nucleus = model.compartments.get_one(id='n')
        mitochondrion = model.compartments.get_one(id='m')
        peroxisome = model.compartments.get_one(id='x')
        er = model.compartments.get_one(id='r')

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

        # Generate response function for the tRNA(s) of each codon and for each amino acid
        trna_kb = cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType,
                                        type=wc_kb.eukaryote.TranscriptType.tRna)
        trna_grouping = {'c': {}, 'm': {}}
        for trna in trna_kb:
            anticodon_prop = trna.properties.get_one(property='anticodon:amino_acid').get_value().split(':')
            codons = ANTICODON_CODON_RECOGNITION_RULES[anticodon_prop[0]]
            for codon in codons:
                if trna.species[0].compartment.id == 'm':
                    if codon in trna_grouping['m']: 
                        trna_grouping['m'][codon]['trna'].append(trna.id)
                        trna_grouping['m'][codon]['anticodon'].append(anticodon_prop[0])
                    else:
                        trna_grouping['m'][codon] = {
                            'trna': [trna.id], 
                            'anticodon': [anticodon_prop[0]], 
                            'aa': anticodon_prop[1], 
                            }
                else:
                    if codon in trna_grouping['c']: 
                        trna_grouping['c'][codon]['trna'].append(trna.id)
                        trna_grouping['c'][codon]['anticodon'].append(anticodon_prop[0])
                    else:
                        trna_grouping['c'][codon] = {
                            'trna': [trna.id],
                            'anticodon': [anticodon_prop[0]], 
                            'aa': anticodon_prop[1], 
                            }

        trna_functions = {'c': {}, 'm': {}}
        for comp, all_trnas in trna_grouping.items():
            for codon, trnas in all_trnas.items():
                compartment = mitochondrion if comp=='m' else cytosol

                # Import cytosolic tRNAs if mitochondrial tRNAs are not detected
                if comp=='m' and all(model.distribution_init_concentrations.get_one(
                  id='dist-init-conc-{}[m]'.format(i)).mean==0 for i in trnas['trna']):
                    trnas['trna'] += trna_grouping['c'][codon]['trna']
                    self._import_cytosolic_trna_into_mitochondria(trna_grouping['c'][codon]['trna'])

                factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                    model, beta, 'translation_{}'.format(compartment.id), 'translation_{}'.format(compartment.id), 
                    compartment, [trnas['trna']])

                objects = {
                        wc_lang.Species: all_species,
                        wc_lang.Parameter: all_parameters,
                        wc_lang.Observable: all_observables,
                        wc_lang.Function: all_volumes,            
                        }

                anticodons = '_'.join(sorted(set(trnas['anticodon'])))
                trna_factor_function = model.functions.get_one(     
                    id='trna_function_{}_{}'.format(anticodons, compartment.id))

                if not trna_factor_function:                   

                    trna_expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                    assert error is None, str(error)

                    trna_factor_function = model.functions.create(     
                        id='trna_function_{}_{}'.format(anticodons, compartment.id),               
                        name='tRNA response function for anticodon(s) {} in {}'.format(
                            anticodons, compartment.name),
                        expression=trna_expression,
                        units=unit_registry.parse_units(''),
                        )
                
                trna_functions[comp][codon] = {
                    'function': trna_factor_function,
                    'aa': trnas['aa'],
                    'objects':objects,
                    }

        aa_functions = {'c': {}, 'm': {}}
        for aa, aa_id in amino_acid_id_conversion.items():
            if aa!='U':
                for compartment in [cytosol, mitochondrion]:
                    factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                        model, beta, 'translation_{}'.format(compartment.id), 'translation_{}'.format(compartment.id), 
                        compartment, [[aa_id]])

                    objects = {
                        wc_lang.Species: all_species,
                        wc_lang.Parameter: all_parameters,
                        wc_lang.Observable: all_observables,
                        wc_lang.Function: all_volumes,            
                        }                                
                    
                    aa_expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                    assert error is None, str(error)

                    aa_functions[compartment.id][aa_id] = {
                        'function': model.functions.create(     
                            id='aminoacid_function_{}_{}'.format(aa_id, compartment.id),               
                            name='response function for amino acid {} in {}'.format(aa_id, compartment.name),
                            expression=aa_expression,
                            units=unit_registry.parse_units(''),
                            ),
                        'objects': objects
                        }
                            
        # Generate response function for each translation initiation factor group
        init_factor_functions = {'c': {}, 'm': {}}
        for comp, factors in {cytosol: cytoplasmic_initiation_factors, mitochondrion: mitochondrial_initiation_factors}.items():
            n = 1
            for factor in factors:
                factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                    model, beta, 'translation_init_{}'.format(comp.id), 'translation_init_{}'.format(comp.id), comp, [factor])

                objects = {
                        wc_lang.Species: all_species,
                        wc_lang.Parameter: all_parameters,
                        wc_lang.Observable: all_observables,
                        wc_lang.Function: all_volumes,            
                        }                                
                    
                expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                assert error is None, str(error)

                init_factor_functions[comp.id][','.join(factor)] = {
                    'function': model.functions.create(     
                        id='translation_init_factor_function_{}_{}'.format(comp.id, n),               
                        name='response function for translation initiation factor {} in {}'.format(n, comp.name),
                        expression=expression,
                        units=unit_registry.parse_units(''),
                        ),
                    'objects': objects}
                n += 1

        # Generate response function for each translation elongation factor group
        el_factor_functions = {'c': {}, 'm': {}}
        for comp, factors in {cytosol: cytoplasmic_elongation_factors, mitochondrion: mitochondrial_elongation_factors}.items():
            n = 1
            for factor in factors:
                factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                    model, beta, 'translation_el_{}'.format(comp.id), 'translation_el_{}'.format(comp.id), comp, [factor])

                objects = {
                        wc_lang.Species: all_species,
                        wc_lang.Parameter: all_parameters,
                        wc_lang.Observable: all_observables,
                        wc_lang.Function: all_volumes,            
                        }                                
                    
                expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                assert error is None, str(error)

                el_factor_functions[comp.id][','.join(factor)] = {
                    'function': model.functions.create(     
                        id='translation_el_factor_function_{}_{}'.format(comp.id, n),               
                        name='response function for translation elongation factor {} in {}'.format(n, comp.name),
                        expression=expression,
                        units=unit_registry.parse_units(''),
                        ),
                    'objects': objects}
                n += 1
        
        # Generate response function for each translocation factor/chaperone group
        trans_factor_functions = {'c': {}, 'm': {}, 'r': {}}
        for comp, factors in {cytosol: cytoplasmic_chaperones, mitochondrion: mitochondrial_chaperones, er: er_chaperones}.items():
            n = 1
            for factor in factors:
                factor_exp, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                    model, beta, 'translocation_{}'.format(comp.id), 'translocation_{}'.format(comp.id), comp, [factor])

                objects = {
                        wc_lang.Species: all_species,
                        wc_lang.Parameter: all_parameters,
                        wc_lang.Observable: all_observables,
                        wc_lang.Function: all_volumes,            
                        }                                
                    
                expression, error = wc_lang.FunctionExpression.deserialize(factor_exp[0], objects)
                assert error is None, str(error)

                trans_factor_functions[comp.id][','.join(factor)] = {
                    'function': model.functions.create(     
                        id='translocation_factor_function_{}_{}'.format(comp.id, n),               
                        name='response function for translocation factor {} in {}'.format(n, comp.name),
                        expression=expression,
                        units=unit_registry.parse_units(''),
                        ),
                    'objects': objects}
                n += 1                
        
        rate_law_no = 0  
        mrna_kbs = [i for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType) \
            if i.type==wc_kb.eukaryote.TranscriptType.mRna]      
        for mrna_kb in mrna_kbs:
            
            mrna_kb_compartment_id = mrna_kb.species[0].compartment.id
            if mrna_kb_compartment_id == 'c':
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
                        
            specific_binding_constant = model.parameters.create(
                id='{}_ribosome_binding_constant'.format(mrna_kb.id),
                type=None,
                units=unit_registry.parse_units('molecule^-2 s^-1'),
                )

            objects = {
                wc_lang.Species: {},
                wc_lang.Parameter: {},
                wc_lang.Observable: {},
                wc_lang.Function: {},            
                }
            expression_terms = []
            for factor in initiation_factors:
                factor_details = init_factor_functions[translation_compartment.id][','.join(factor)]
                expression_terms.append(factor_details['function'].id)
                for cl, dictionary in objects.items():
                    dictionary.update(factor_details['objects'][cl])
                objects[wc_lang.Function][factor_details['function'].id] = factor_details['function']

            start_codon = gvar.protein_aa_usage[mrna_kb.protein.id]['start_codon'].replace('U', 'T')
            start_aa_met = amino_acid_id_conversion[gvar.protein_aa_usage[mrna_kb.protein.id]['start_aa']]
            matched_trnas = [trna_functions[translation_compartment.id][start_codon]]            
            for codon_info in matched_trnas:
                expression_terms.append(codon_info['function'].id)
                objects[wc_lang.Function][codon_info['function'].id] = codon_info['function']                
                for cl, dictionary in objects.items():
                    dictionary.update(codon_info['objects'][cl])                 
                           
            expression_terms.append(aa_functions[translation_compartment.id][
                start_aa_met]['function'].id)
            objects[wc_lang.Function][aa_functions[translation_compartment.id][
                start_aa_met]['function'].id] = aa_functions[translation_compartment.id][
                start_aa_met]['function']
            for cl, dictionary in objects.items():
                dictionary.update(aa_functions[translation_compartment.id][
                    start_aa_met]['objects'][cl])    
            
            objects[wc_lang.Species][ribosome_complex.id] = ribosome_complex
            objects[wc_lang.Species][self._allowable_queue_len[mrna_kb.id][0].id]= self._allowable_queue_len[mrna_kb.id][0]

            objects[wc_lang.Parameter][max_bool.id] = max_bool
            objects[wc_lang.Parameter][min_bool.id] = min_bool
            objects[wc_lang.Parameter][specific_binding_constant.id] = specific_binding_constant
            
            expression = '{} * {} * max(min({} , {}) , {}) * {} * 2**{}'.format(
                specific_binding_constant.id,
                ribosome_complex.id,
                self._allowable_queue_len[mrna_kb.id][0].id,
                max_bool.id,
                min_bool.id,
                ' * '.join(expression_terms),
                len(expression_terms),
                )

            init_rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, objects)
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

            objects = {
                wc_lang.Species: {},
                wc_lang.Parameter: {},
                wc_lang.Observable: {},
                wc_lang.Function: {},            
                }
            expression_terms = []
            for factor in elongation_factors:
                factor_details = el_factor_functions[translation_compartment.id][','.join(factor)]
                expression_terms.append(factor_details['function'].id)
                for cl, dictionary in objects.items():
                    dictionary.update(factor_details['objects'][cl])
                objects[wc_lang.Function][factor_details['function'].id] = factor_details['function']
            
            if codon_table == 1:
                codon_id = 1
            else:
                codon_id = codon_table[mrna_kb.protein.id]
            coding_rna_seq, _, _ = mrna_kb.protein.get_seq_and_start_codon(table=codon_id, cds=cds)
            codon_seq = str(coding_rna_seq).upper().replace('U','T')
            all_codons = sorted(set([codon_seq[i * 3:(i + 1) * 3] for i in range((len(codon_seq) + 3 - 1) // 3 )][1:]))
            for i in all_codons:
                if len(i)==3:
                    if i in UNCONDITIONAL_STOP_CODON:
                        pass
                    elif i in CONDITIONAL_STOP_CODON and mrna_kb.gene.id not in selenoproteome:
                        pass    
                    else:
                        if translation_compartment.id == 'c':    
                            matched_trnas = [trna_functions[translation_compartment.id][i]]
                        else:
                            if i in trna_functions[translation_compartment.id]:
                                matched_trnas = [trna_functions[translation_compartment.id][i]]
                            else:
                                cytosolic_trna_ids = trna_grouping['c'][i]['trna']
                                self._import_cytosolic_trna_into_mitochondria(cytosolic_trna_ids)

                                factor_exp, all_species, all_parameters, all_volumes, all_observables = \
                                    utils.gen_response_functions(model, beta, 'translation_m', 'translation_m', 
                                    translation_compartment, [cytosolic_trna_ids])

                                added_objects = {
                                    wc_lang.Species: all_species,
                                    wc_lang.Parameter: all_parameters,
                                    wc_lang.Observable: all_observables,
                                    wc_lang.Function: all_volumes,            
                                    }

                                anticodons = '_'.join(sorted(set(trna_grouping['c'][i]['anticodon'])))
                                trna_factor_function = model.functions.get_one(     
                                    id='trna_function_{}_m'.format(anticodons))

                                if not trna_factor_function:

                                    trna_expression, error = wc_lang.FunctionExpression.deserialize(
                                        factor_exp[0], added_objects)
                                    assert error is None, str(error)
                                
                                    trna_factor_function = model.functions.create(     
                                        id='trna_function_{}_m'.format(anticodons),               
                                        name='tRNA response function for anticodon(s) {} in {}'.format(
                                            anticodons, translation_compartment.name),
                                        expression=trna_expression,
                                        units=unit_registry.parse_units(''),
                                        )
                                
                                trna_functions['m'][i] = {
                                    'function': trna_factor_function,
                                    'aa': trna_functions['c'][i]['aa'],
                                    'objects': added_objects,
                                    }
                                matched_trnas = [trna_functions[translation_compartment.id][i]]       
                            
                        for codon_info in matched_trnas:    
                            expression_terms.append(codon_info['function'].id)
                            objects[wc_lang.Function][codon_info['function'].id] = codon_info['function']                        
                            for cl, dictionary in objects.items():
                                dictionary.update(codon_info['objects'][cl])                    

            selcys = 0
            for key, value in gvar.protein_aa_usage[mrna_kb.protein.id].items():
                aa_id = 'S' if key=='U' else key                
                if aa_id in amino_acid_id_conversion and value:
                    selcys += 1 if key=='U' else 0
                    aa_met = amino_acid_id_conversion[aa_id]
                    if aa_met==start_aa_met and value-1==0: 
                        pass
                    else:    
                        if aa_functions[translation_compartment.id][aa_met][
                            'function'].id not in expression_terms:
                            
                            expression_terms.append(aa_functions[translation_compartment.id][
                                aa_met]['function'].id)
                            objects[wc_lang.Function][aa_functions[translation_compartment.id][
                                aa_met]['function'].id] = aa_functions[translation_compartment.id][
                                aa_met]['function']
                            for cl, dictionary in objects.items():
                                dictionary.update(aa_functions[translation_compartment.id][
                                    aa_met]['objects'][cl])
            
            other_mets = [['gtp'], ['atp']] + ([['selnp']] if selcys else [])
            expressions, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                model, beta, elongation_reaction.id, 'translation_elongation', translation_compartment, other_mets)
            expression_terms += expressions
            objects[wc_lang.Species].update(all_species)
            objects[wc_lang.Parameter].update(all_parameters)
            objects[wc_lang.Function].update(all_volumes)
            objects[wc_lang.Observable].update(all_observables)
            
            ribo_bound_species = model.species_types.get_one(id='ribo_bound_{}'.format(mrna_kb.id)).species[0]
            objects[wc_lang.Species][ribo_bound_species.id] = ribo_bound_species

            k_cat_elongation = model.parameters.create(
                id='k_cat_{}'.format(elongation_reaction.id),
                type=wc_ontology['WC:k_cat'],
                units=unit_registry.parse_units('molecule^-1 s^-1'),
                )
            objects[wc_lang.Parameter][k_cat_elongation.id] = k_cat_elongation

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
                
                if target_compartment_id=='n':
                    energy_compartment = nucleus
                    energy_reactant = 'gtp'
                    chaperones = []
                elif target_compartment_id=='m':
                    energy_compartment = mitochondrion
                    energy_reactant = 'atp'
                    chaperones = mitochondrial_chaperones
                    chaperone_compartment = mitochondrion                    
                elif target_compartment_id=='x':
                    energy_compartment = peroxisome
                    energy_reactant = 'atp'
                    chaperones = []
                else:
                    energy_compartment = cytosol
                    energy_reactant = 'gtp'
                    chaperones = er_chaperones
                    chaperone_compartment = er

                objects = {
                    wc_lang.Species: {},
                    wc_lang.Parameter: {},
                    wc_lang.Observable: {},
                    wc_lang.Function: {},            
                    }
                expression_terms = []
                for factor in chaperones:
                    factor_details = trans_factor_functions[chaperone_compartment.id][','.join(factor)]
                    expression_terms.append(factor_details['function'].id)
                    for cl, dictionary in objects.items():
                        dictionary.update(factor_details['objects'][cl])
                    objects[wc_lang.Function][factor_details['function'].id] = factor_details['function']    

                expressions, all_species, all_parameters, all_volumes, all_observables = utils.gen_response_functions(
                    model, beta, reaction.id, 'translocation', energy_compartment, [[energy_reactant]])
                expression_terms += expressions
                objects[wc_lang.Species].update(all_species)
                objects[wc_lang.Parameter].update(all_parameters)
                objects[wc_lang.Function].update(all_volumes)
                objects[wc_lang.Observable].update(all_observables)

                k_cat_translocation = model.parameters.create(
                    id='k_cat_{}'.format(reaction.id),
                    type=wc_ontology['WC:k_cat'],
                    units=unit_registry.parse_units('molecule^-1 s^-1'),
                    )
                objects[wc_lang.Parameter][k_cat_translocation.id] = k_cat_translocation

                protein_species = model.species_types.get_one(
                    id=reaction.id.split('_')[1]).species.get_one(compartment=translation_compartment)
                objects[wc_lang.Species][protein_species.id] = protein_species

                volume = translation_compartment.init_density.function_expressions[0].function
                objects[wc_lang.Function][volume.id] = volume

                expression = '{} * {} * {} * 2**{}'.format(
                    k_cat_translocation.id,
                    protein_species.id,
                    ' * '.join(expression_terms),
                    len(expression_terms),
                    )

                trans_rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, objects)
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
        peroxisome = model.compartments.get_one(id='x')
        er = model.compartments.get_one(id='r')

        init_compartment_volumes = {
            nucleus.id: nucleus.init_volume.mean * nucleus.init_density.value,
            mitochondrion.id: mitochondrion.init_volume.mean * mitochondrion.init_density.value,
            cytosol.id: cytosol.init_volume.mean * cytosol.init_density.value,
            peroxisome.id: peroxisome.init_volume.mean * peroxisome.init_density.value,
            er.id: er.init_volume.mean * er.init_density.value,
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
            
            if mrna_kb_compartment_id == 'c':   

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
        determined_init_kcat = []
        undetermined_init_kcat = []
        determined_el_kcat = []
        undetermined_el_kcat = []
        determined_transloc_kcat = []
        undetermined_transloc_kcat = []
        for mrna_kb in mrna_kbs:

            mrna_kb_compartment_id = mrna_kb.species[0].compartment.id
            ribosome_complex = cytoplasmic_ribosome_species if mrna_kb_compartment_id == 'c' else mitochondrial_ribosome_species

            protein_model = model.species_types.get_one(id=mrna_kb.protein.id)
            
            complex_model_stoic = {model.species_types.get_one(id=i.id):j.coefficient for i in cell.species_types.get(
                __type=wc_kb.core.ComplexSpeciesType) for j in i.subunits if j.species_type==mrna_kb.protein}

            total_concentration = sum([i.distribution_init_concentration.mean for i in protein_model.species]) + \
                sum([i.distribution_init_concentration.mean*v for k,v in complex_model_stoic.items() for i in k.species \
                    if i.distribution_init_concentration])
            half_life = mrna_kb.protein.properties.get_one(property='half-life').get_value()

            average_rate = utils.calc_avg_syn_rate(
                total_concentration, half_life, mean_doubling_time)

            # Calibrate initiation reaction
            init_reaction = model.reactions.get_one(id='translation_initiation_' + mrna_kb.id)

            init_species_counts = {}
                        
            init_species_counts[ribosome_complex.id] = ribosome_complex.distribution_init_concentration.mean
            init_species_counts[self._allowable_queue_len[mrna_kb.id][0].id] = self._allowable_queue_len[mrna_kb.id][1]

            for species in init_reaction.rate_laws[0].expression.species:
                init_species_counts[species.id] = species.distribution_init_concentration.mean
                model_Km = model.parameters.get_one(
                        id='K_m_{}_{}'.format(init_reaction.id, species.species_type.id))
                if model_Km:
                    if species.distribution_init_concentration.mean:                    
                        model_Km.value = beta * species.distribution_init_concentration.mean \
                            / Avogadro.value / species.compartment.init_volume.mean
                        model_Km.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                            beta, species.species_type.id, species.compartment.name)
                    else:
                        model_Km.value = 1e-05
                        model_Km.comments = 'The value was assigned to 1e-05 because the concentration of ' +\
                            '{} in {} was zero'.format(species.species_type.id, species.compartment.name)    

            for func in init_reaction.rate_laws[0].expression.functions:
                for species in func.expression.species:
                    init_species_counts[species.id] = species.distribution_init_concentration.mean
                for obs in func.expression.observables:
                    for species in obs.expression.species:    
                        init_species_counts[species.id] = species.distribution_init_concentration.mean

            model_kcat = model.parameters.get_one(id='{}_ribosome_binding_constant'.format(mrna_kb.id))
            
            if average_rate:            
                model_kcat.value = 1.
                eval_rate_law = init_reaction.rate_laws[0].expression._parsed_expression.eval({
                    wc_lang.Species: init_species_counts,
                    wc_lang.Compartment: init_compartment_volumes,
                })
                if eval_rate_law:
                    model_kcat.value = average_rate / eval_rate_law
                    determined_init_kcat.append(model_kcat.value)
                else:
                    undetermined_init_kcat.append(model_kcat)    
            else:          
                model_kcat.value = 0.

            # Calibrate elongation reaction
            el_reaction = model.reactions.get_one(id='translation_elongation_' + mrna_kb.id)
            
            el_species_counts = {}
            
            for species in el_reaction.rate_laws[0].expression.species:
                el_species_counts[species.id] = species.distribution_init_concentration.mean
                model_Km = model.parameters.get_one(
                        id='K_m_{}_{}'.format(el_reaction.id, species.species_type.id))
                if model_Km:
                    if species.distribution_init_concentration.mean:                     
                        model_Km.value = beta * species.distribution_init_concentration.mean \
                            / Avogadro.value / species.compartment.init_volume.mean
                        model_Km.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                            beta, species.species_type.id, species.compartment.name)
                    else:
                        model_Km.value = 1e-05
                        model_Km.comments = 'The value was assigned to 1e-05 because the concentration of ' +\
                        '{} in {} was zero'.format(species.species_type.id, species.compartment.name)    
                    
            for func in el_reaction.rate_laws[0].expression.functions:                
                for species in func.expression.species:
                    el_species_counts[species.id] = species.distribution_init_concentration.mean
                for obs in func.expression.observables:
                    for species in obs.expression.species:    
                        el_species_counts[species.id] = species.distribution_init_concentration.mean

            model_kcat = model.parameters.get_one(id='k_cat_{}'.format(el_reaction.id))

            if average_rate:            
                model_kcat.value = 1.
                eval_rate_law = el_reaction.rate_laws[0].expression._parsed_expression.eval({
                    wc_lang.Species: el_species_counts,
                    wc_lang.Compartment: init_compartment_volumes,
                    })
                if eval_rate_law:
                    model_kcat.value = average_rate / eval_rate_law
                    determined_el_kcat.append(model_kcat.value)
                else:
                    undetermined_el_kcat.append(model_kcat)     
            else:          
                model_kcat.value = 0.

            # Calibrate translocation reaction
            conc_per_comp = {}
            for protein in protein_model.species:
                conc_per_comp[protein.compartment] = protein.distribution_init_concentration.mean
            for cplx_st, stoic in complex_model_stoic.items():
                for cplx_species in cplx_st.species:
                    if cplx_species.distribution_init_concentration:
                        if cplx_species.compartment in conc_per_comp:
                            conc_per_comp[cplx_species.compartment] += stoic * \
                                cplx_species.distribution_init_concentration.mean
                        else:        
                            conc_per_comp[cplx_species.compartment] = stoic * \
                                cplx_species.distribution_init_concentration.mean
            
            translation_compartment = cytosol if mrna_kb_compartment_id == 'c' else mitochondrion

            for compartment, trans_reaction in self._translocation_reactions[mrna_kb].items():
                trans_species_counts = {}
                
                for species in trans_reaction.rate_laws[0].expression.species:
                    trans_species_counts[species.id] = species.distribution_init_concentration.mean
                    model_Km = model.parameters.get_one(
                            id='K_m_{}_{}'.format(trans_reaction.id, species.species_type.id))
                    if model_Km:
                        if species.distribution_init_concentration.mean:                        
                            model_Km.value = beta * species.distribution_init_concentration.mean \
                                / Avogadro.value / species.compartment.init_volume.mean
                            model_Km.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                                beta, species.species_type.id, species.compartment.name)
                        else:
                            model_Km.value = 1e-05
                            model_Km.comments = 'The value was assigned to 1e-05 because the concentration of ' +\
                                '{} in {} was zero'.format(species.species_type.id, species.compartment.name)    
                            
                for func in trans_reaction.rate_laws[0].expression.functions:
                    for species in func.expression.species:
                        trans_species_counts[species.id] = species.distribution_init_concentration.mean
                    for obs in func.expression.observables:
                        for species in obs.expression.species:    
                            trans_species_counts[species.id] = species.distribution_init_concentration.mean

                model_kcat = model.parameters.get_one(id='k_cat_{}'.format(trans_reaction.id))

                if average_rate:            
                    model_kcat.value = 1.
                    eval_rate_law = trans_reaction.rate_laws[0].expression._parsed_expression.eval({
                        wc_lang.Species: trans_species_counts,
                        wc_lang.Compartment: init_compartment_volumes,
                        })
                    if eval_rate_law:
                        model_kcat.value = conc_per_comp[compartment] / conc_per_comp[translation_compartment] * \
                            average_rate / eval_rate_law
                        determined_transloc_kcat.append(model_kcat.value)
                    else:
                        undetermined_transloc_kcat.append(model_kcat)         
                else:          
                    model_kcat.value = 0.            
            
        median_init_kcat = numpy.median(determined_init_kcat)
        if not numpy.isnan(median_init_kcat): 
            for model_kcat in undetermined_init_kcat:
                model_kcat.value = median_init_kcat
                model_kcat.comments = 'Set to the median value because it could not be determined from data'
        else:
            for model_kcat in undetermined_init_kcat:
                model_kcat.value = 1.
                model_kcat.comments = 'Set to 1 because it could not be determined from median value'                 

        median_el_kcat = numpy.median(determined_el_kcat)
        if not numpy.isnan(median_el_kcat):
            for model_kcat in undetermined_el_kcat:
                model_kcat.value = median_el_kcat
                model_kcat.comments = 'Set to the median value because it could not be determined from data'
        else:
            for model_kcat in undetermined_el_kcat:
                model_kcat.value = 1.
                model_kcat.comments = 'Set to 1 because it could not be determined from median value'       

        median_transloc_kcat = numpy.median(determined_transloc_kcat)
        if not numpy.isnan(median_transloc_kcat):
            for model_kcat in undetermined_transloc_kcat:
                model_kcat.value = median_transloc_kcat
                model_kcat.comments = 'Set to the median value because it could not be determined from data'
        else:
            for model_kcat in undetermined_transloc_kcat:
                model_kcat.value = 1.
                model_kcat.comments = 'Set to 1 because it could not be determined from median value'

    def _import_cytosolic_trna_into_mitochondria(self, cytosolic_trna_ids):
        """ Create reactions and rate laws for importing cytosolic tRNAs into the mitochondria.
            The concentrations of imported tRNAs in the mitochondria are set based on
            the provided fraction and the rates of transport are calibrated accordingly to achieve
            steady-states

            Args:
                cytosolic_trna_ids (:obj:`list`): list of species type IDs of cytosolic tRNAa to be
                    imported into the mitochondria
        """
        kb = self.knowledge_base            
        model = self.model
        submodel = self.submodel
        beta = self.options['beta']
        mitochondrial_cytosolic_trna_partition = self.options['mitochondrial_cytosolic_trna_partition']
        mitochondrial_exosome = self.options['mitochondrial_exosome']

        cytosol = model.compartments.get_one(id='c')
        mitochondria = model.compartments.get_one(id='m')
        
        rna_deg_submodel = model.submodels.get_one(id='rna_degradation')
        exosome_species = model.species_types.get_one(
            name=mitochondrial_exosome).species.get_one(compartment=mitochondria)

        Avogadro = model.parameters.get_or_create(
            id='Avogadro',
            type=None,
            value=scipy.constants.Avogadro,
            units=unit_registry.parse_units('molecule mol^-1'))       

        for trna_id in cytosolic_trna_ids:

            trna_species_type = model.species_types.get_one(id=trna_id)

            if not trna_species_type.species.get_one(compartment=mitochondria):

                cyto_trna_species = trna_species_type.species.get_one(compartment=cytosol)
                total_conc = cyto_trna_species.distribution_init_concentration.mean
                cyto_trna_species.distribution_init_concentration.mean = total_conc * \
                    (1 - mitochondrial_cytosolic_trna_partition)
                cyto_trna_species.distribution_init_concentration.comments = \
                    'Value is adjusted to account for import into the mitochondria'

                mito_trna_species = model.species.get_or_create(
                    species_type=trna_species_type, compartment=mitochondria)
                mito_trna_species.id = mito_trna_species.gen_id()

                conc_model = model.distribution_init_concentrations.create(
                    species=mito_trna_species,
                    mean=total_conc * mitochondrial_cytosolic_trna_partition,
                    units=unit_registry.parse_units('molecule'),
                    comments='Value is set to {} of the total cellular concentration'.format(
                        mitochondrial_cytosolic_trna_partition),
                    )
                conc_model.id = conc_model.gen_id()

                # Generate import reaction
                import_reaction = model.reactions.create(
                    submodel=self.submodel, id='trna_import_{}'.format(trna_id),
                    name='import of {} into the mitochondria'.format(trna_id),
                    reversible=False)
                import_reaction.participants.append(
                    cyto_trna_species.species_coefficients.get_or_create(coefficient=-1))
                import_reaction.participants.append(
                    mito_trna_species.species_coefficients.get_or_create(coefficient=1))

                # Generate rate law for import reaction
                import_constant = model.parameters.create(
                    id='{}_import_constant'.format(trna_id),
                    type=None,
                    units=unit_registry.parse_units('molecule^-1 s^-1'),
                    )
                expression = '{} * {}'.format(import_constant.id, cyto_trna_species.id)

                import_rate_law_expression, error = wc_lang.RateLawExpression.deserialize(expression, {
                    wc_lang.Species: {cyto_trna_species.id: cyto_trna_species},
                    wc_lang.Parameter: {import_constant.id: import_constant},
                    })
                assert error is None, str(error)

                import_rate_law = model.rate_laws.create(
                    direction=wc_lang.RateLawDirection.forward,
                    type=None,
                    expression=import_rate_law_expression,
                    reaction=import_reaction,
                    units=unit_registry.parse_units('s^-1'),                
                    )
                import_rate_law.id = import_rate_law.gen_id()

                # Calibrate import rate
                rna_kb = kb.cell.species_types.get_one(id=trna_id)
                half_life = rna_kb.properties.get_one(property='half-life').get_value()
                mean_doubling_time = model.parameters.get_one(id='mean_doubling_time').value
                average_rate = utils.calc_avg_syn_rate(
                    total_conc, half_life, mean_doubling_time)
                
                import_constant.value = 1.
                eval_rate_law = import_rate_law_expression._parsed_expression.eval({
                    wc_lang.Species: {cyto_trna_species.id: total_conc * \
                        (1 - mitochondrial_cytosolic_trna_partition)}
                    })
                if eval_rate_law:
                    import_constant.value = mitochondrial_cytosolic_trna_partition / \
                        (1 - mitochondrial_cytosolic_trna_partition) * \
                        average_rate / eval_rate_law
                else:        
                    import_constant.value = 0.

                # Generate degradation reaction for imported tRNA
                metabolic_participants = ['amp', 'cmp', 'gmp', 'ump', 'h2o', 'h']
                metabolites = {}
                for met in metabolic_participants:
                    met_species_type = model.species_types.get_one(id=met)
                    metabolites[met] = met_species_type.species.get_or_create(
                            compartment=mitochondria, model=model)

                reaction = model.reactions.get_or_create(submodel=rna_deg_submodel, 
                    id='degradation_{}_{}'.format(trna_id, mitochondria.id))
                reaction.name = 'degradation of {} in mitochondria'.format(trna_species_type.name)
                
                if trna_id in gvar.transcript_ntp_usage:
                    ntp_count = gvar.transcript_ntp_usage[trna_id]
                else:
                    seq = rna_kb.get_seq()
                    ntp_count = gvar.transcript_ntp_usage[trna_id] = {
                        'A': seq.upper().count('A'),
                        'C': seq.upper().count('C'),
                        'G': seq.upper().count('G'),
                        'U': seq.upper().count('U'),
                        'len': len(seq)
                        }                
                # Adding participants to LHS
                reaction.participants.append(mito_trna_species.species_coefficients.get_or_create(
                    coefficient=-1))
                reaction.participants.append(metabolites['h2o'].species_coefficients.get_or_create(
                    coefficient=-(ntp_count['len']-1)))
                # Adding participants to RHS
                reaction.participants.append(metabolites['amp'].species_coefficients.get_or_create(
                    coefficient=ntp_count['A']))
                reaction.participants.append(metabolites['cmp'].species_coefficients.get_or_create(
                    coefficient=ntp_count['C']))
                reaction.participants.append(metabolites['gmp'].species_coefficients.get_or_create(
                    coefficient=ntp_count['G']))
                reaction.participants.append(metabolites['ump'].species_coefficients.get_or_create(
                    coefficient=ntp_count['U']))
                reaction.participants.append(metabolites['h'].species_coefficients.get_or_create(
                    coefficient=ntp_count['len']-1))            
                
                # Generate rate law for degradation reaction of imported tRNA
                rate_law_exp, _ = utils.gen_michaelis_menten_like_rate_law(
                    model, reaction, modifiers=[exosome_species], 
                    exclude_substrates=[metabolites['h2o']])
                
                rate_law = model.rate_laws.create(
                    direction=wc_lang.RateLawDirection.forward,
                    type=None,
                    expression=rate_law_exp,
                    reaction=reaction,
                    )
                rate_law.id = rate_law.gen_id()
                
                # Calibrate degradation reaction of imported tRNA
                cyto_deg_reaction = rna_deg_submodel.reactions.get_one(id='degradation_' + rna_kb.id)
                cyto_species_counts = {species.id: species.distribution_init_concentration.mean \
                    for species in cyto_deg_reaction.rate_laws[0].expression.species}
                cyto_deg_rate = cyto_deg_reaction.rate_laws[0].expression._parsed_expression.eval({
                    wc_lang.Species: cyto_species_counts,
                    wc_lang.Compartment: {
                        cytosol.id: cytosol.init_volume.mean * \
                            cytosol.init_density.value}
                    })                
                total_deg_rate = utils.calc_avg_deg_rate(total_conc, half_life)
                mito_deg_rate = total_deg_rate - cyto_deg_rate

                mito_species_counts = {exosome_species.id: exosome_species.distribution_init_concentration.mean}
                for species in reaction.get_reactants():

                    mito_species_counts[species.id] = species.distribution_init_concentration.mean

                    if model.parameters.get(id='K_m_{}_{}'.format(reaction.id, species.species_type.id)):
                        model_Km = model.parameters.get_one(
                            id='K_m_{}_{}'.format(reaction.id, species.species_type.id))
                        if species.distribution_init_concentration.mean:
                            model_Km.value = beta * species.distribution_init_concentration.mean \
                                / Avogadro.value / species.compartment.init_volume.mean
                            model_Km.comments = 'The value was assumed to be {} times the concentration of {} in {}'.format(
                                beta, species.species_type.id, species.compartment.name)
                        else:
                            model_Km.value = 1e-05
                            model_Km.comments = 'The value was assigned to 1e-05 because the concentration of ' +\
                                '{} in {} was zero'.format(species.species_type.id, species.compartment.name)

                model_kcat = model.parameters.get_one(id='k_cat_{}'.format(reaction.id))

                if mito_deg_rate:            
                    model_kcat.value = 1.
                    eval_rate_law = reaction.rate_laws[0].expression._parsed_expression.eval({
                        wc_lang.Species: mito_species_counts,
                        wc_lang.Compartment: {
                            mitochondria.id: mitochondria.init_volume.mean * \
                                mitochondria.init_density.value}
                        })
                    if eval_rate_law:
                        model_kcat.value = mito_deg_rate / eval_rate_law
                    else:
                        model_kcat.value = 0.  
                else:          
                    model_kcat.value = 0.
