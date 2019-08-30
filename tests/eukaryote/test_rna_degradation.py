""" Tests of rna degradation submodel generation
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-06-11
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import rna_degradation
from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import math
import os
import scipy.constants
import shutil
import tempfile
import unittest
import wc_lang
import wc_kb
import wc_kb_gen


class RnaDegradationSubmodelGeneratorTestCase(unittest.TestCase):

    def setUp(self):

        # Create KB content
        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>chr1\nTTTATGACTCTAGTTTAT\n'
                    '>chrM\nTTTATGACTC TAGTTTAT\n')

        self.kb = wc_kb.KnowledgeBase()
        cell = self.kb.cell = wc_kb.Cell()

        nucleus = cell.compartments.create(id='n')
        mito = cell.compartments.create(id='m')

        chr1 = wc_kb.core.DnaSpeciesType(cell=cell, id='chr1', sequence_path=self.sequence_path)
        gene1 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene1', polymer=chr1, start=1, end=18)
        exon1 = wc_kb.eukaryote.GenericLocus(start=4, end=18)
        transcript1 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans1', 
            name='transcript1', gene=gene1, exons=[exon1])
        transcript1_half_life = wc_kb.core.SpeciesTypeProperty(property='half_life', species_type=transcript1, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        transcript1_spec = wc_kb.core.Species(species_type=transcript1, compartment=nucleus)
        transcript1_conc = wc_kb.core.Concentration(cell=cell, species=transcript1_spec, value=10.)

        chrM = wc_kb.core.DnaSpeciesType(cell=cell, id='chrM', sequence_path=self.sequence_path)
        gene2 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene2', polymer=chrM, start=1, end=18)
        exon2 = wc_kb.eukaryote.GenericLocus(start=1, end=10)
        transcript2 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans2', 
            name='transcript2', gene=gene2, exons=[exon2])
        transcript2_half_life = wc_kb.core.SpeciesTypeProperty(property='half_life', species_type=transcript2, 
            value='15000.0', value_type=wc_ontology['WC:float'])
        transcript2_spec = wc_kb.core.Species(species_type=transcript2, compartment=mito)
        transcript2_conc = wc_kb.core.Concentration(cell=cell, species=transcript2_spec, value=10.)

        transcript3 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans3', 
            name='transcript3', gene=gene2, exons=[exon2])
        transcript3_half_life = wc_kb.core.SpeciesTypeProperty(property='half_life', species_type=transcript3, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        transcript3_spec = wc_kb.core.Species(species_type=transcript3, compartment=mito)
        transcript3_conc = wc_kb.core.Concentration(cell=cell, species=transcript3_spec, value=10.)                   

        # Create initial model content
        self.model = model = wc_lang.Model()
        
        model.parameters.create(id='Avogadro', value = scipy.constants.Avogadro,
                                units = unit_registry.parse_units('molecule mol^-1'))

        compartments = {'n': ('nucleus', 5E-14), 'm': ('mitochondria', 2.5E-14), 'c': ('cytoplasm', 9E-14)}
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

        for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType):
            model_species_type = model.species_types.create(id=i.id)
            model_compartment = model.compartments.get_one(id='m' if 'M' in i.gene.polymer.id else 'n')
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=model_species, 
                mean=10., units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()

        complexes = {'complex1': ('Exosome', ['c', 'n']), 'complex2': ('Exosome variant', ['c', 'n']), 'complex3': ('Mitochondrial Exosome', ['m']),
            'complex4': ('Mitochondrial Exosome variant', ['m'])}
        for k, v in complexes.items():
            model_species_type = model.species_types.get_or_create(id=k, name=v[0])
            for comp in v[1]:
                model_compartment = model.compartments.get_one(id=comp)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, 
                    mean=100., units=unit_registry.parse_units('molecule'))
                conc_model.id = conc_model.gen_id()

        metabolic_participants = ['amp', 'cmp', 'gmp', 'ump', 'h2o', 'h']
        for i in metabolic_participants:
            model_species_type = model.species_types.create(id=i)
            for c in ['n', 'm', 'c']:
                model_compartment = model.compartments.get_one(id=c)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, 
                    mean=1500., units=unit_registry.parse_units('molecule'))
                conc_model.id = conc_model.gen_id()

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)            

    def test_methods(self):

        gen = rna_degradation.RnaDegradationSubmodelGenerator(self.kb, self.model, options={
            'rna_exo_pair': {'trans1': 'Exosome', 'trans2': 'Mitochondrial Exosome', 
            'trans3': 'Mitochondrial Exosome'}
            })
        gen.run()

        # Test gen_reactions
        self.assertEqual([i.id for i in self.model.submodels], ['rna_degradation'])
        self.assertEqual(sorted([i.id for i in self.model.reactions]), 
            sorted(['degradation_trans1', 'degradation_trans2', 'degradation_trans3']))
        self.assertEqual(sorted([i.name for i in self.model.reactions]), 
            sorted(['degradation of transcript1', 'degradation of transcript2', 'degradation of transcript3']))
        self.assertEqual(set([i.submodel.id for i in self.model.reactions]), set(['rna_degradation']))
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='degradation_trans1').participants}, 
            {'amp[c]': 4, 'cmp[c]': 2, 'gmp[c]': 2, 'ump[c]': 7, 'h[c]': 14, 'h2o[c]': -14, 'trans1[n]': -1})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='degradation_trans2').participants}, 
            {'amp[m]': 2, 'cmp[m]': 2, 'gmp[m]': 1, 'ump[m]': 5, 'h[m]': 9, 'h2o[m]': -9, 'trans2[m]': -1})
        
        # Test gen_rate_laws
        self.assertEqual(len(self.model.rate_laws), 3)
        self.assertEqual(self.model.rate_laws.get_one(id='degradation_trans1-forward').expression.expression,
            'k_cat_degradation_trans1 * complex1[c] * '
            '(trans1[n] / (trans1[n] + K_m_degradation_trans1_trans1 * Avogadro * volume_n))')
        self.assertEqual(self.model.rate_laws.get_one(id='degradation_trans2-forward').expression.expression,
            'k_cat_degradation_trans2 * complex3[m] * '
            '(trans2[m] / (trans2[m] + K_m_degradation_trans2_trans2 * Avogadro * volume_m))')
        self.assertEqual(self.model.rate_laws.get_one(id='degradation_trans3-forward').expression.expression,
            'k_cat_degradation_trans3 * complex3[m] * '
            '(trans3[m] / (trans3[m] + K_m_degradation_trans3_trans3 * Avogadro * volume_m))')
        
        # Test calibrate_submodel
        self.assertEqual(self.model.parameters.get_one(id='K_m_degradation_trans2_trans2').value, 10/scipy.constants.Avogadro/2.5E-14)
        self.assertEqual(self.model.parameters.get_one(id='k_cat_degradation_trans1').value, math.log(2)/36000*10/(0.5*100))
        self.assertEqual(self.model.parameters.get_one(id='k_cat_degradation_trans2').value, math.log(2)/15000*10/(0.5*100))
