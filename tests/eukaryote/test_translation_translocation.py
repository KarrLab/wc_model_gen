""" Tests of translation-translocation submodel generation
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-12-09
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import translation_translocation
from wc_onto import onto as wc_ontology
from wc_utils.util.chem import EmpiricalFormula
from wc_utils.util.units import unit_registry
import wc_model_gen.global_vars as gvar
import math
import os
import scipy.constants
import shutil
import tempfile
import unittest
import wc_lang
import wc_kb
import wc_kb_gen

class TranslationTranslocationSubmodelGeneratorTestCase(unittest.TestCase):

    def setUp(self):

    	# Create KB content
        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>chr1\nATGGCGTGCGATGAT\n'
                    '>chrM\nNGCGTGCATGGATGAT\n')

        self.kb = wc_kb.KnowledgeBase()
        cell = self.kb.cell = wc_kb.Cell()

        nucleus = cell.compartments.create(id='n')
        cytoplasm = cell.compartments.create(id='c')
        mito = cell.compartments.create(id='m')
        
        chr1 = wc_kb.core.DnaSpeciesType(cell=cell, id='chr1', sequence_path=self.sequence_path)
        gene1 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene1', polymer=chr1, start=1, end=12)
        chrM = wc_kb.core.DnaSpeciesType(cell=cell, id='chrM', sequence_path=self.sequence_path)
        geneM = wc_kb.eukaryote.GeneLocus(cell=cell, id='geneM', polymer=chrM, start=2, end=16)

        locus1 = wc_kb.eukaryote.GenericLocus(start=1, end=9)
        transcript1 = wc_kb.eukaryote.TranscriptSpeciesType(id='trans1', cell=cell, gene=gene1, exons=[locus1], type=wc_kb.eukaryote.TranscriptType.mRna)
        transcript1_spec = wc_kb.core.Species(species_type=transcript1, compartment=nucleus)
        transcript1_conc = wc_kb.core.Concentration(cell=cell, species=transcript1_spec, value=10.)
        prot1 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot1', name='protein1', transcript=transcript1, coding_regions=[locus1])
        prot1_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot1, 
            value='40000.0', value_type=wc_ontology['WC:float'])

        locus2 = wc_kb.eukaryote.GenericLocus(start=1, end=9)
        transcript2 = wc_kb.eukaryote.TranscriptSpeciesType(id='trans2', cell=cell, gene=gene1, exons=[locus2], type=wc_kb.eukaryote.TranscriptType.mRna)
        transcript2_spec = wc_kb.core.Species(species_type=transcript2, compartment=nucleus)
        transcript2_conc = wc_kb.core.Concentration(cell=cell, species=transcript2_spec, value=10.)
        prot2 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot2', name='protein2', transcript=transcript2, coding_regions=[locus2])
        prot2_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot2, 
            value='10000.0', value_type=wc_ontology['WC:float'])

        locus3 = wc_kb.eukaryote.GenericLocus(start=4, end=9)
        transcript3 = wc_kb.eukaryote.TranscriptSpeciesType(id='trans3', cell=cell, gene=gene1, exons=[locus3], type=wc_kb.eukaryote.TranscriptType.rRna)
        
        locusM = wc_kb.eukaryote.GenericLocus(start=8, end=16)
        transcriptM = wc_kb.eukaryote.TranscriptSpeciesType(id='transM', cell=cell, gene=geneM, exons=[locusM], type=wc_kb.eukaryote.TranscriptType.mRna)
        transcriptM_spec = wc_kb.core.Species(species_type=transcriptM, compartment=mito)
        transcriptM_conc = wc_kb.core.Concentration(cell=cell, species=transcriptM_spec, value=10.)
        protM = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='protM', name='proteinM', transcript=transcriptM, coding_regions=[locusM])
        protM_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=protM, 
            value='25000.0', value_type=wc_ontology['WC:float'])

        complexM = wc_kb.core.ComplexSpeciesType(id='complexM', cell=cell, subunits=[wc_kb.core.SpeciesTypeCoefficient(species_type=protM,
            coefficient=2)])

        trnas = {'mt_trnaM': ('M', 'CAT', 'm'), 'mt_trnaA': ('A', 'CGC', 'm'), 'mt_trnaC': ('C', 'GCA', 'm'), 'mt_trnaD': ('D', 'ATC', 'm'), 
            'trnaM': ('M', 'CAT', 'c'), 'trnaA1': ('A', 'CGC', 'c'), 'trnaA2': ('A', 'CGC', 'c'), 'trnaC': ('C', 'GCA', 'c'), 'trnaD': ('D', 'ATC', 'c')}
        for trna_id, (aa_id, anticodon, comp) in trnas.items():
            trna_species_type = wc_kb.eukaryote.TranscriptSpeciesType(id=trna_id, cell=cell, type=wc_kb.eukaryote.TranscriptType.tRna)
            trna_compartment = cytoplasm if comp=='c' else mito
            trna_spec = wc_kb.core.Species(species_type=trna_species_type, compartment=trna_compartment)            
            trna_property = wc_kb.core.SpeciesTypeProperty(property='anticodon:amino_acid', species_type=trna_species_type, 
                value='{}:{}'.format(anticodon, aa_id), value_type=wc_ontology['WC:string'])

        # Create initial model content
        self.model = model = wc_lang.Model()

        model.parameters.create(id='mean_doubling_time', value=20*3600, units=unit_registry.parse_units('s'))

        compartments = {'n': ('nucleus', 5E-14), 'c': ('cytoplasm', 1E-13), 'm': ('mitochondria', 2.5E-14), 
            'x': ('peroxisome', 5E-15), 'r': ('endoplasmic reticulum', 1.5E-14), 'c_m': ('membrane', 5E-16)}
        for k, v in compartments.items():
            init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], 
                    mean=v[1], std=0)
            c = model.compartments.create(id=k, name=v[0], init_volume=init_volume)
            c.init_density = model.parameters.create(id='density_' + k, value=1000, 
                units=unit_registry.parse_units('g l^-1'))
            volume = model.functions.create(id='volume_' + k, units=unit_registry.parse_units('l'))
            volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
            assert error is None, str(error)

        transcripts = {'trans1': ('c', [5, 15]), 'trans2': ('c', [5, 15]), 'transM': ('m', [3, 9])}
        for Id, details in transcripts.items():
            model_species_type = model.species_types.create(id=Id)
            model_compartment = model.compartments.get_one(id=details[0])
            model_species = model.species.create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=model_species, 
                mean=details[1][0], units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()
            site_species_type = model.species_types.create(id='{}_ribosome_binding_site'.format(Id))
            site_compartment = model.compartments.get_one(id=details[0])
            site_species = model.species.create(species_type=site_species_type, compartment=site_compartment)
            site_species.id = site_species.gen_id()
            site_conc = model.distribution_init_concentrations.create(species=site_species,
                mean=details[1][1], units=unit_registry.parse_units('molecule'))
            site_conc.id = site_conc.gen_id()            

        proteins = {'prot1': ['n', 'c', 'x', 'r', 'm', 'c_m'], 'prot2': ['c'], 'protM': ['m']}
        for k, v in proteins.items():
            kb_protein = cell.species_types.get_one(id=k)
            model_species_type = model.species_types.create(id=kb_protein.id, name=kb_protein.name, type=wc_ontology['WC:protein'])
            for comp in v:
                model_compartment = model.compartments.get_one(id=comp)
                model_species = model.species.create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, 
                    mean=10, units=unit_registry.parse_units('molecule'))
                conc_model.id = conc_model.gen_id()
        model.distribution_init_concentrations.get_one(id='dist-init-conc-prot2[c]').mean = 0.
        model.distribution_init_concentrations.get_one(id='dist-init-conc-prot1[c_m]').mean = 5        
        
        complexes = {'comp_1': ('ribosome', ['m', 'c']), 'comp_2': ('init_factor', ['m', 'c']), 'comp_3': ('el_factor1', ['m', 'c']), 
            'comp_4': ('el_factor2', ['m', 'c']), 'comp_5': ('chaperone', ['n', 'c', 'x', 'r', 'm']), 
            'complexM': ('complexM', ['m'])}
        for k, v in complexes.items():
            model_species_type = model.species_types.create(id=k, name=v[0], type=wc_ontology['WC:pseudo_species'],
            	structure = wc_lang.ChemicalStructure(
                empirical_formula = EmpiricalFormula('H'),
                molecular_weight = 1.018,
                charge = 1))
            for comp in v[1]:
                model_compartment = model.compartments.get_one(id=comp)
                model_species = model.species.create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, 
                    mean=5, units=unit_registry.parse_units('molecule'))
                conc_model.id = conc_model.gen_id()

        metabolic_participants = ['Met', 'Ala', 'Cys', 'Asp', 'atp', 'adp', 'amp', 'gtp','gdp','pi', 'h2o', 'h']
        for i in metabolic_participants:
            model_species_type = model.species_types.create(id=i, type=wc_ontology['WC:metabolite'],
                structure = wc_lang.ChemicalStructure(
                empirical_formula = EmpiricalFormula('H'),
                molecular_weight = 1.018,
                charge = 1))
            for c in compartments:                            
                model_compartment = model.compartments.get_one(id=c)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, 
                    mean=15, units=unit_registry.parse_units('molecule'))
                conc_model.id = conc_model.gen_id()                    

        trna_model = {'mt_trnaM': 'm', 'mt_trnaA':'m', 'mt_trnaC': 'm', 'mt_trnaD': 'm', 
            'trnaM': 'c','trnaA1': 'c', 'trnaA2': 'c', 'trnaC': 'c', 'trnaD': 'c'}
        for Id, c in trna_model.items():
            model_species_type = model.species_types.create(id=Id)
            model_compartment = model.compartments.get_one(id=c)
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()  
            conc_model = model.distribution_init_concentrations.create(species=model_species, 
                mean=2, units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)
        gvar.protein_aa_usage = {}
        gvar.transcript_ntp_usage = {}          

    def test_methods(self):           

        gen = translation_translocation.TranslationTranslocationSubmodelGenerator(self.kb, self.model, options={
            'cytoplasmic_ribosome': 'ribosome',
            'mitochondrial_ribosome': 'ribosome',
            'cytoplasmic_initiation_factors': [['init_factor']],
            'mitochondrial_initiation_factors': [['init_factor']],
            'cytoplasmic_elongation_factors': [['prot2', 'el_factor2']],
            'mitochondrial_elongation_factors': [['el_factor1'], ['protM']],
            'cytoplasmic_chaperones': [['chaperone']],
            'mitochondrial_chaperones': [['chaperone']],
            'er_chaperones': [['chaperone']],
            'amino_acid_id_conversion': {
                'M': 'Met',
                'A': 'Ala',
                'C': 'Cys',
                'D': 'Asp',
                },
            'cds': False,
            'polysome_fraction': {'trans1': 0.4, 'trans2': 0.2, 'transM': 0.4},
            })
        gen.run()

        model = self.model

        self.assertEqual(gvar.transcript_ntp_usage['trans1'], {'A': 1, 'U': 2, 'G': 4, 'C': 2, 'len': 9})
        self.assertEqual(gvar.protein_aa_usage['prot1'], {'A': 1, 'C': 1, 'D': 0, 'M': 1, 'len': 3, '*': 0})

        # Test gen_reactions
        self.assertEqual([i.id for i in self.model.submodels], ['translation_translocation'])

        # initiation
        self.assertEqual(model.species_types.get_one(id='ribo_bound_trans1').structure.empirical_formula, EmpiricalFormula('H2'))
        self.assertEqual(model.species_types.get_one(id='ribo_bound_trans1').structure.molecular_weight, 2.036)
        self.assertEqual(model.species_types.get_one(id='ribo_bound_trans1').structure.charge, 2)        
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translation_initiation_trans1').participants},
            {'trans1_ribosome_binding_site[c]': -1, 'comp_1[c]': -1, 'Met[c]': -1, 'h2o[c]': -5, 'atp[c]': -2, 'gtp[c]': -2,
             'ribo_bound_trans1[c]': 1, 'h[c]': 5, 'amp[c]': 1, 'adp[c]': 1, 'gdp[c]': 2, 'pi[c]': 5})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translation_initiation_transM').participants},
            {'transM_ribosome_binding_site[m]': -1, 'comp_1[m]': -1, 'Met[m]': -1, 'h2o[m]': -5, 'atp[m]': -2, 'gtp[m]': -2,
             'ribo_bound_transM[m]': 1, 'h[m]': 5, 'amp[m]': 1, 'adp[m]': 1, 'gdp[m]': 2, 'pi[m]': 5})

        # elongation        
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translation_elongation_trans1').participants}, 
            {'ribo_bound_trans1[c]': -1, 'gtp[c]': -3, 'atp[c]': -2, 'h2o[c]': -5, 'Ala[c]': -1, 'Cys[c]': -1, 
            'comp_1[c]': 1, 'trans1_ribosome_binding_site[c]': 1, 'prot1[c]': 1, 'amp[c]': 2, 'gdp[c]': 3, 'pi[c]': 7, 'h[c]': 3})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translation_elongation_transM').participants}, 
            {'ribo_bound_transM[m]': -1, 'gtp[m]': -3, 'atp[m]': -2, 'h2o[m]': -5, 'Asp[m]': -2, 
            'comp_1[m]': 1, 'transM_ribosome_binding_site[m]': 1, 'protM[m]': 1, 'amp[m]': 2, 'gdp[m]': 3, 'pi[m]': 7, 'h[m]': 3})
        
        # translocation
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translocation_prot1_c_to_n').participants}, 
            {'prot1[c]': -1, 'gtp[n]': -1, 'h2o[n]': -1, 'prot1[n]': 1, 'gdp[n]': 1, 'pi[n]': 1, 'h[n]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translocation_prot1_c_to_x').participants}, 
            {'prot1[c]': -1, 'atp[x]': -1, 'h2o[x]': -1, 'prot1[x]': 1, 'adp[x]': 1, 'pi[x]': 1, 'h[x]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translocation_prot1_c_to_r').participants}, 
            {'prot1[c]': -1, 'gtp[r]': -1, 'h2o[r]': -1, 'prot1[r]': 1, 'gdp[r]': 1, 'pi[r]': 1, 'h[r]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translocation_prot1_c_to_m').participants}, 
            {'prot1[c]': -1, 'atp[m]': -1, 'h2o[m]': -1, 'prot1[m]': 1, 'adp[m]': 1, 'pi[m]': 1, 'h[m]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translocation_prot1_c_to_c_m').participants}, 
            {'prot1[c]': -1, 'gtp[r]': -1, 'h2o[r]': -1, 'prot1[c_m]': 1, 'gdp[r]': 1, 'pi[r]': 1, 'h[r]': 1})
        self.assertEqual(model.reactions.get_one(id='translocation_protM_m_to_m'), None)

        # Test gen_rate_laws
        self.assertEqual(len(model.rate_laws), 11)
        self.assertEqual(len(model.observables), 2)
        self.assertEqual(len(model.functions), 30) # 6 volume + 8 trna + 8 aa + 2 init + 3 elo + 3 trans
        self.assertEqual(model.observables.get_one(id='translation_c_factors_c_1').expression.expression, 'trnaA1[c] + trnaA2[c]')
        self.assertEqual(model.observables.get_one(id='translation_el_c_factors_c_1').expression.expression, 'prot2[c] + comp_4[c]')
        self.assertEqual(model.functions.get_one(id='trna_function_GCG_c').expression.expression, 
        	'(translation_c_factors_c_1 / (translation_c_factors_c_1 + K_m_translation_c_translation_c_factors_c_1 * Avogadro * volume_c))')
        self.assertEqual(model.functions.get_one(id='trna_function_GCG_m').expression.expression, 
        	'(mt_trnaA[m] / (mt_trnaA[m] + K_m_translation_m_mt_trnaA * Avogadro * volume_m))')
        self.assertEqual(model.functions.get_one(id='aminoacid_function_Ala_c').expression.expression, 
        	'(Ala[c] / (Ala[c] + K_m_translation_c_Ala * Avogadro * volume_c))')
        self.assertEqual(model.functions.get_one(id='aminoacid_function_Asp_m').expression.expression, 
        	'(Asp[m] / (Asp[m] + K_m_translation_m_Asp * Avogadro * volume_m))')
        self.assertEqual(model.functions.get_one(id='translation_init_factor_function_c_1').expression.expression, 
        	'(comp_2[c] / (comp_2[c] + K_m_translation_init_c_comp_2 * Avogadro * volume_c))')
        self.assertEqual(model.functions.get_one(id='translation_init_factor_function_m_1').expression.expression, 
        	'(comp_2[m] / (comp_2[m] + K_m_translation_init_m_comp_2 * Avogadro * volume_m))')
        self.assertEqual(model.functions.get_one(id='translation_el_factor_function_c_1').expression.expression, 
        	'(translation_el_c_factors_c_1 / (translation_el_c_factors_c_1 + K_m_translation_el_c_translation_el_c_factors_c_1 * Avogadro * volume_c))')
        self.assertEqual(model.functions.get_one(id='translation_el_factor_function_m_1').expression.expression, 
        	'(comp_3[m] / (comp_3[m] + K_m_translation_el_m_comp_3 * Avogadro * volume_m))')
        self.assertEqual(model.functions.get_one(id='translation_el_factor_function_m_2').expression.expression, 
        	'(protM[m] / (protM[m] + K_m_translation_el_m_protM * Avogadro * volume_m))')
        self.assertEqual(model.functions.get_one(id='translocation_factor_function_c_1').expression.expression, 
        	'(comp_5[c] / (comp_5[c] + K_m_translocation_c_comp_5 * Avogadro * volume_c))')
        self.assertEqual(model.functions.get_one(id='translocation_factor_function_m_1').expression.expression, 
        	'(comp_5[m] / (comp_5[m] + K_m_translocation_m_comp_5 * Avogadro * volume_m))')
        self.assertEqual(model.functions.get_one(id='translocation_factor_function_r_1').expression.expression, 
        	'(comp_5[r] / (comp_5[r] + K_m_translocation_r_comp_5 * Avogadro * volume_r))')

        # Initiation
        self.assertEqual(model.rate_laws.get_one(id='translation_initiation_trans1-forward').expression.expression,
            'trans1_ribosome_binding_constant * comp_1[c] * '
            'max(min(trans1_ribosome_binding_site[c] , max_bool_substance) , min_bool_substance) * '
            'translation_init_factor_function_c_1 * '
            'trna_function_AUG_c * '
            'aminoacid_function_Met_c * ' 
            '2**3')
        self.assertEqual(model.rate_laws.get_one(id='translation_initiation_transM-forward').expression.expression,
            'transM_ribosome_binding_constant * comp_1[m] * '
            'max(min(transM_ribosome_binding_site[m] , max_bool_substance) , min_bool_substance) * '
            'translation_init_factor_function_m_1 * '
            'trna_function_AUG_m * '
            'aminoacid_function_Met_m * ' 
            '2**3')

        # Elongation
        self.assertEqual(model.rate_laws.get_one(id='translation_elongation_trans1-forward').expression.expression,
            'k_cat_translation_elongation_trans1 * ribo_bound_trans1[c] * '
            'translation_el_factor_function_c_1 * '
            'trna_function_GCG_c * '
            'aminoacid_function_Ala_c * '
            'trna_function_UGC_c * '
            'aminoacid_function_Cys_c * '
            '(gtp[c] / (gtp[c] + K_m_translation_elongation_trans1_gtp * Avogadro * volume_c)) * '
            '(atp[c] / (atp[c] + K_m_translation_elongation_trans1_atp * Avogadro * volume_c)) * '
            '2**7')
        self.assertEqual(model.rate_laws.get_one(id='translation_elongation_transM-forward').expression.expression,
            'k_cat_translation_elongation_transM * ribo_bound_transM[m] * '
            'translation_el_factor_function_m_1 * '
            'translation_el_factor_function_m_2 * '
            'trna_function_GAU_m * '
            'aminoacid_function_Asp_m * '
            '(gtp[m] / (gtp[m] + K_m_translation_elongation_transM_gtp * Avogadro * volume_m)) * '
            '(atp[m] / (atp[m] + K_m_translation_elongation_transM_atp * Avogadro * volume_m)) * '
            '2**6')

        # Translocation
        self.assertEqual(model.rate_laws.get_one(id='translocation_prot1_c_to_m-forward').expression.expression,
            'k_cat_translocation_prot1_c_to_m * prot1[c] * '
            'translocation_factor_function_m_1 * '
            '(atp[m] / (atp[m] + K_m_translocation_prot1_c_to_m_atp * Avogadro * volume_m)) * '
            '2**2')
        self.assertEqual(model.rate_laws.get_one(id='translocation_prot1_c_to_c_m-forward').expression.expression,
            'k_cat_translocation_prot1_c_to_c_m * prot1[c] * '
            'translocation_factor_function_r_1 * '
            '(gtp[r] / (gtp[r] + K_m_translocation_prot1_c_to_c_m_gtp * Avogadro * volume_r)) * '
            '2**2')

        for law in model.rate_laws:
            self.assertEqual(law.validate(), None)

        # Test calibrate_submodel
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-ribo_bound_trans1[c]').mean, 2)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-ribo_bound_transM[m]').mean, 2)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-comp_1[c]').mean, 2)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-comp_1[m]').mean, 3)

        # Initiation
        self.assertEqual(model.parameters.get_one(id='K_m_translation_init_c_comp_2').value, 5/scipy.constants.Avogadro/1E-13)
        self.assertEqual(model.parameters.get_one(id='K_m_translation_init_c_comp_2').comments, 
            'The value was assumed to be 1.0 times the concentration of comp_2 in cytoplasm')
        self.assertEqual(model.parameters.get_one(id='trans1_ribosome_binding_constant').value, math.log(2)*(1/(20*3600) + 1/40000)*55/2)
        self.assertEqual(model.parameters.get_one(id='trans2_ribosome_binding_constant').value, 0.)
        
        # Elongation
        self.assertEqual(model.parameters.get_one(id='K_m_translation_el_c_translation_el_c_factors_c_1').value, (0+5)/scipy.constants.Avogadro/1E-13)
        self.assertEqual(model.parameters.get_one(id='K_m_translation_el_c_translation_el_c_factors_c_1').comments, 
        	'The value was assumed to be 1.0 times the value of translation_el_c_factors_c_1')
        self.assertEqual(model.parameters.get_one(id='k_cat_translation_elongation_transM').value, math.log(2)*(1/(20*3600) + 1/25000)*20/2)
        self.assertEqual(model.parameters.get_one(id='k_cat_translation_elongation_trans2').value, 0.)

        # Translocation
        self.assertEqual(model.parameters.get_one(id='K_m_translocation_m_comp_5').value, 5/scipy.constants.Avogadro/2.5E-14)
        self.assertEqual(model.parameters.get_one(id='K_m_translocation_m_comp_5').comments, 
        	'The value was assumed to be 1.0 times the concentration of comp_5 in mitochondria')
        self.assertEqual(model.parameters.get_one(id='K_m_translocation_prot1_c_to_c_m_gtp').value, 15/scipy.constants.Avogadro/1.5E-14)
        self.assertEqual(model.parameters.get_one(id='k_cat_translocation_prot1_c_to_m').value, math.log(2)*(1/(20*3600) + 1/40000)*55/10)
        self.assertEqual(model.parameters.get_one(id='k_cat_translocation_prot1_c_to_c_m').value, math.log(2)*(1/(20*3600) + 1/40000)*55/10*5/10)

    def test_global_vars(self):
        gvar.protein_aa_usage = {'prot1': {'M': 2, 'A': 4, 'C': 2, 'D': 1, 'len': 7, '*': 1}}
        gen = translation_translocation.TranslationTranslocationSubmodelGenerator(self.kb, self.model, options={
            'cytoplasmic_ribosome': 'ribosome',
            'mitochondrial_ribosome': 'ribosome',
            'cytoplasmic_initiation_factors': [['init_factor']],
            'mitochondrial_initiation_factors': [['init_factor']],
            'cytoplasmic_elongation_factors': [['el_factor1', 'el_factor2']],
            'mitochondrial_elongation_factors': [['el_factor1'], ['el_factor2']],
            'cytoplasmic_chaperones': [['chaperone']],
            'mitochondrial_chaperones': [['chaperone']],
            'er_chaperones': [['chaperone']],
            'amino_acid_id_conversion': {
                'M': 'Met',
                'A': 'Ala',
                'C': 'Cys',
                'D': 'Asp',
                },
            'cds': False,
            'polysome_fraction': {'trans1': 0.4, 'trans2': 0.2, 'transM': 0.4},
            })
        gen.run()

        model = self.model

        self.assertEqual(gvar.protein_aa_usage['prot1'], {'M':2, 'A': 4, 'C': 2, 'D': 1, 'len': 7, '*': 1})

        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translation_elongation_trans1').participants}, 
            {'ribo_bound_trans1[c]': -1, 'gtp[c]': -7, 'atp[c]': -6, 'h2o[c]': -13, 'Ala[c]': -4, 'Cys[c]': -2, 'Asp[c]': -1, 'Met[c]': -1,
            'comp_1[c]': 1, 'trans1_ribosome_binding_site[c]': 1, 'prot1[c]': 1, 'amp[c]': 6, 'gdp[c]': 7, 'pi[c]': 19, 'h[c]': 7})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='translation_elongation_transM').participants}, 
            {'ribo_bound_transM[m]': -1, 'gtp[m]': -3, 'atp[m]': -2, 'h2o[m]': -5, 'Asp[m]': -2, 
            'comp_1[m]': 1, 'transM_ribosome_binding_site[m]': 1, 'protM[m]': 1, 'amp[m]': 2, 'gdp[m]': 3, 'pi[m]': 7, 'h[m]': 3})
