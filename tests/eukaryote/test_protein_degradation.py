""" Tests of protein degradation submodel generation
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-06-11
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import protein_degradation
from wc_onto import onto as wc_ontology
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


class ProteinDegradationSubmodelGeneratorTestCase(unittest.TestCase):

    def setUp(self):

        # Create KB content
        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>chr1\nGCGTGCGATGAT\n'
                    '>chrM\nNGCGTGCGATGAT\n')

        self.kb = wc_kb.KnowledgeBase()
        cell = self.kb.cell = wc_kb.Cell()

        chr1 = wc_kb.core.DnaSpeciesType(cell=cell, id='chr1', sequence_path=self.sequence_path)
        gene1 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene1', polymer=chr1, start=1, end=12)
        chrM = wc_kb.core.DnaSpeciesType(cell=cell, id='chrM', sequence_path=self.sequence_path)
        geneM = wc_kb.eukaryote.GeneLocus(cell=cell, id='geneM', polymer=chrM, start=2, end=13)

        locus1 = wc_kb.eukaryote.GenericLocus(start=1, end=6)
        transcript1 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, gene=gene1, exons=[locus1])
        prot1 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot1', name='protein1', transcript=transcript1, coding_regions=[locus1])
        prot1_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot1, 
            value='40000.0', value_type=wc_ontology['WC:float'])
        
        locus2 = wc_kb.eukaryote.GenericLocus(start=4, end=9)
        transcript2 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, gene=gene1, exons=[locus2])
        prot2 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot2', name='protein2', transcript=transcript2, coding_regions=[locus2])
        prot2_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot2, 
            value='40000.0', value_type=wc_ontology['WC:float'])
                
        locus3 = wc_kb.eukaryote.GenericLocus(start=7, end=12)
        transcript3 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, gene=gene1, exons=[locus3])
        prot3 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot3', name='protein3', transcript=transcript3, coding_regions=[locus3])
        prot3_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot3, 
            value='25000.0', value_type=wc_ontology['WC:float'])
        
        locusM = wc_kb.eukaryote.GenericLocus(start=8, end=13)
        transcriptM = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, gene=geneM, exons=[locusM])
        protM = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='protM', name='proteinM', transcript=transcriptM, coding_regions=[locusM])
        protM_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=protM, 
            value='25000.0', value_type=wc_ontology['WC:float'])
                
        # Create initial model content
        self.model = model = wc_lang.Model()
        
        model.parameters.create(id='Avogadro', value = scipy.constants.Avogadro,
                                units = unit_registry.parse_units('molecule mol^-1'))

        compartments = {'n': ('nucleus', 5E-14), 'c': ('cytoplasm', 1E-13),
            'm': ('mitochondria', 2.5E-14), 'l': ('lysosome', 1.5E-14), 'c_m': ('membrane', 5E-15)}
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

        proteins = {'prot1': ['n', 'c'], 'prot2': ['c', 'l'], 'prot3': ['l'], 'protM': ['m', 'c_m']}
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
        model.distribution_init_concentrations.get_one(id='dist-init-conc-protM[m]').mean = 0.        

        complexes = {'comp_1': ('26S proteasome', ['n', 'c']), 'comp_2': ('lonp1', ['m']), 'comp_3': ('clpp', ['m']), 
            'comp_4': ('cathepsin B', ['l']), 'comp_5': ('cathepsin D', ['l'])}
        for k, v in complexes.items():
            model_species_type = model.species_types.create(id=k, name=v[0], type=wc_ontology['WC:pseudo_species'])
            for comp in v[1]:
                model_compartment = model.compartments.get_one(id=comp)
                model_species = model.species.create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, 
                    mean=2, units=unit_registry.parse_units('molecule'))
                conc_model.id = conc_model.gen_id()

        metabolic_participants = ['Ala', 'Cys', 'Asp', 'Glu', 'h2o']
        metabolic_compartments = ['l', 'm']
        for i in metabolic_participants:
            for c in metabolic_compartments:
                model_species_type = model.species_types.get_or_create(id=i, type=wc_ontology['WC:metabolite'])            
                model_compartment = model.compartments.get_one(id=c)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()                    

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)
        gvar.protein_aa_usage = {}            

    def test_methods(self):

        gen = protein_degradation.ProteinDegradationSubmodelGenerator(self.kb, self.model, options={
            'compartment_proteasomes': {
                'n': ['26S proteasome'], 
                'c': ['prot2'], 
                'l': ['cathepsin B', 'cathepsin D'], 
                'm': ['lonp1', 'protM'],
                },
            'amino_acid_id_conversion': {
                'A': 'Ala',
                'C': 'Cys',
                'D': 'Asp',
                },
            'cds': False,
            })
        gen.run()

        self.assertEqual(gvar.protein_aa_usage['prot1'], {'A': 1, 'C': 1, 'D': 0, 'len': 2, '*': 0})

        # Test gen_reactions
        self.assertEqual([i.id for i in self.model.submodels], ['protein_degradation'])
        self.assertEqual(sorted([i.id for i in self.model.reactions]), 
            sorted(['prot1_n_degradation', 'prot1_c_degradation', 'prot2_c_degradation', 
                'prot2_l_degradation', 'prot3_l_degradation', 'protM_m_degradation', 'protM_c_m_degradation']))
        self.assertEqual(sorted([i.name for i in self.model.reactions]), 
            sorted(['Degradation of prot1 of nucleus', 'Degradation of prot1 of cytoplasm', 'Degradation of prot2 of cytoplasm', 
                'Degradation of prot2 of lysosome', 'Degradation of prot3 of lysosome', 'Degradation of protM of mitochondria',
                'Degradation of protM of membrane']))
        self.assertEqual(set([i.submodel.id for i in self.model.reactions]), set(['protein_degradation']))
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='prot1_n_degradation').participants}, 
            {'prot1[n]': -1, 'h2o[n]': -1, 'Ala[n]': 1, 'Cys[n]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='prot1_c_degradation').participants}, 
            {'prot1[c]': -1, 'h2o[c]': -1, 'Ala[c]': 1, 'Cys[c]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='prot2_c_degradation').participants}, 
            {'prot2[c]': -1, 'h2o[c]': -1, 'Asp[c]': 1, 'Cys[c]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='prot2_l_degradation').participants}, 
            {'prot2[l]': -1, 'h2o[l]': -1, 'Asp[l]': 1, 'Cys[l]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='prot3_l_degradation').participants}, 
            {'prot3[l]': -1, 'h2o[l]': -1, 'Asp[l]': 2})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='protM_m_degradation').participants}, 
            {'protM[m]': -1, 'h2o[m]': -1, 'Asp[m]': 2})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='protM_c_m_degradation').participants}, 
            {'protM[c_m]': -1, 'h2o[l]': -1, 'Asp[l]': 2})
        
        # Test gen_rate_laws
        self.assertEqual(len(self.model.rate_laws), 7)
        self.assertEqual(len(self.model.observables), 2)
        self.assertEqual(self.model.observables.get_one(id='total_proteasomes_l').expression.expression, 'comp_4[l] + comp_5[l]')
        self.assertEqual(self.model.observables.get_one(id='total_proteasomes_m').expression.expression, 'comp_2[m] + protM[m]')
        self.assertEqual(self.model.rate_laws.get_one(id='prot1_n_degradation-forward').expression.expression,
            'k_cat_prot1_n_degradation * comp_1[n] * '
            '(prot1[n] / (prot1[n] + K_m_prot1_n_degradation_prot1 * Avogadro * volume_n))')
        self.assertEqual(self.model.rate_laws.get_one(id='protM_m_degradation-forward').expression.expression,
            'k_cat_protM_m_degradation * total_proteasomes_m * '
            '(protM[m] / (protM[m] + K_m_protM_m_degradation_protM * Avogadro * volume_m))')
        self.assertEqual(self.model.rate_laws.get_one(id='protM_c_m_degradation-forward').expression.expression,
            'k_cat_protM_c_m_degradation * total_proteasomes_l * '
            '(protM[c_m] / (protM[c_m] + K_m_protM_c_m_degradation_protM * Avogadro * volume_c_m))')

        for law in self.model.rate_laws:
            self.assertEqual(law.validate(), None)
        
        # Test calibrate_submodel
        self.assertEqual(self.model.parameters.get_one(id='K_m_prot1_n_degradation_prot1').value, 10/scipy.constants.Avogadro/5E-14)
        self.assertEqual(self.model.parameters.get_one(id='K_m_prot1_n_degradation_prot1').comments, 
            'The value was assumed to be 1.0 times the concentration of prot1 in nucleus')
        self.assertEqual(self.model.parameters.get_one(id='K_m_protM_c_m_degradation_protM').value, 10/scipy.constants.Avogadro/5E-15)
        self.assertEqual(self.model.parameters.get_one(id='K_m_protM_c_m_degradation_protM').comments, 
            'The value was assumed to be 1.0 times the concentration of protM in membrane before transport to lysosome')
        self.assertEqual(self.model.parameters.get_one(id='k_cat_prot1_n_degradation').value, math.log(2)/40000*10/(0.5*2))
        self.assertEqual(self.model.parameters.get_one(id='k_cat_protM_c_m_degradation').value, math.log(2)/25000*10/(0.5*4))

        self.assertEqual(self.model.parameters.get_one(id='K_m_protM_m_degradation_protM').value, 0.5*(10/scipy.constants.Avogadro/5E-14 + 10/scipy.constants.Avogadro/1.5E-14))
        self.assertEqual(self.model.parameters.get_one(id='K_m_protM_m_degradation_protM').comments, 
            'Set to the median value because protein concentration was zero')
        self.assertEqual(self.model.parameters.get_one(id='k_cat_protM_m_degradation').value, 0.5*(math.log(2)/40000*10/(0.5*4) + math.log(2)/25000*10/(0.5*4)))
        self.assertEqual(self.model.parameters.get_one(id='k_cat_protM_m_degradation').comments, 
            'Set to the median value because it could not be determined from data')

    def test_global_vars(self):
        gvar.protein_aa_usage = {'prot1': {'A': 4, 'C': 2, 'D': 1, 'len': 7, '*': 1}}
        gen = protein_degradation.ProteinDegradationSubmodelGenerator(self.kb, self.model, options={
            'compartment_proteasomes': {
                'n': ['26S proteasome'], 
                'c': ['26S proteasome'], 
                'l': ['cathepsin B', 'cathepsin D'], 
                'm': ['lonp1', 'clpp'],
                },
            'amino_acid_id_conversion': {
                'A': 'Ala',
                'C': 'Cys',
                'D': 'Asp',
                },
            'cds': False,
            })
        gen.run()

        self.assertEqual(gvar.protein_aa_usage['prot1'], {'A': 4, 'C': 2, 'D': 1, 'len': 7, '*': 1})

        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='prot1_n_degradation').participants}, 
            {'prot1[n]': -1, 'h2o[n]': -6, 'Ala[n]': 4, 'Cys[n]': 2, 'Asp[n]':1})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='protM_c_m_degradation').participants}, 
            {'protM[c_m]': -1, 'h2o[l]': -1, 'Asp[l]': 2})
        