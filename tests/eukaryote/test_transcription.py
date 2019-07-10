""" Tests of transcription submodel generation

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-05-21
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import transcription
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

class TranscriptionSubmodelGeneratorTestCase(unittest.TestCase):

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
        gene1 = wc_kb.eukaryote_schema.GeneLocus(cell=cell, id='gene1', polymer=chr1, start=1, end=18)
        exon1 = wc_kb.eukaryote_schema.GenericLocus(start=4, end=18)
        transcript1 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans1', 
            name='transcript1', gene=gene1, exons=[exon1])
        transcript1_half_life = wc_kb.core.SpeciesTypeProperty(property='half_life', species_type=transcript1, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        transcript1_spec = wc_kb.core.Species(species_type=transcript1, compartment=nucleus)
        transcript1_conc = wc_kb.core.Concentration(cell=cell, species=transcript1_spec, value=10.)

        chrM = wc_kb.core.DnaSpeciesType(cell=cell, id='chrM', sequence_path=self.sequence_path)
        gene2 = wc_kb.eukaryote_schema.GeneLocus(cell=cell, id='gene2', polymer=chrM, start=1, end=18)
        exon2 = wc_kb.eukaryote_schema.GenericLocus(start=1, end=10)
        transcript2 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans2', 
            name='transcript2', gene=gene2, exons=[exon2])
        transcript2_half_life = wc_kb.core.SpeciesTypeProperty(property='half_life', species_type=transcript2, 
            value='15000.0', value_type=wc_ontology['WC:float'])
        transcript2_spec = wc_kb.core.Species(species_type=transcript2, compartment=mito)
        transcript2_conc = wc_kb.core.Concentration(cell=cell, species=transcript2_spec, value=10.)

        transcript3 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans3', 
            name='transcript3', gene=gene2, exons=[exon2])
        transcript3_half_life = wc_kb.core.SpeciesTypeProperty(property='half_life', species_type=transcript3, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        transcript3_spec = wc_kb.core.Species(species_type=transcript3, compartment=mito)
        transcript3_conc = wc_kb.core.Concentration(cell=cell, species=transcript3_spec, value=10.)                   

        # Create initial model content
        self.model = model = wc_lang.Model()
        
        model.parameters.create(id='mean_doubling_time', value=20*3600, units=unit_registry.parse_units('s'))
        model.parameters.create(id='Avogadro', value = scipy.constants.Avogadro,
                                units = unit_registry.parse_units('molecule mol^-1'))

        compartments = {'n': ('nucleus', 5E-14), 'm': ('mitochondria', 2.5E-14)}
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

        for i in cell.species_types.get(__type=wc_kb.eukaryote_schema.TranscriptSpeciesType):
            model_species_type = model.species_types.create(id=i.id)
            model_compartment = model.compartments.get_one(id='m' if 'M' in i.gene.polymer.id else 'n')
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=model_species, 
                mean=10., units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()

        complexes = {'complex1': ('RNA Polymerase I','n'), 'complex2': ('RNA Polymerase I', 'n'), 
            'complex3': ('RNA Polymerase II', 'm'), 'complex4': ('RNA Polymerase III', 'm')}
        for k, v in complexes.items():
            model_species_type = model.species_types.create(id=k, name=v[0])
            model_compartment = model.compartments.get_one(id=v[1])
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=model_species, 
                mean=100., units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()

        metabolic_participants = ['atp', 'ctp', 'gtp', 'utp', 'ppi']
        for i in metabolic_participants:
            model_species_type = model.species_types.create(id=i)
            for c in ['n', 'm']:
                model_compartment = model.compartments.get_one(id=c)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, 
                    mean=1500., units=unit_registry.parse_units('molecule'))
                conc_model.id = conc_model.gen_id()

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)            

    def test_methods(self):            

        gen = transcription.TranscriptionSubmodelGenerator(self.kb, self.model, options={
            'rna_pol_pair': {'trans1': ['RNA Polymerase I'], 'trans2': ['RNA Polymerase II'], 
                            'trans3': ['RNA Polymerase II', 'RNA Polymerase III']}
            })
        gen.run()  

        # Test gen_reactions
        self.assertEqual([i.id for i in self.model.submodels], ['transcription'])
        self.assertEqual(sorted([i.id for i in self.model.reactions]), 
            sorted(['transcription_trans1', 'transcription_trans2', 'transcription_trans3']))
        self.assertEqual(sorted([i.name for i in self.model.reactions]), 
            sorted(['transcription of transcript1', 'transcription of transcript2', 'transcription of transcript3']))
        self.assertEqual(set([i.submodel.id for i in self.model.reactions]), set(['transcription']))
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='transcription_trans1').participants}, 
            {'atp[n]': -4, 'ctp[n]': -2, 'gtp[n]': -2, 'utp[n]': -7, 'ppi[n]': 14, 'trans1[n]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='transcription_trans2').participants}, 
            {'atp[m]': -2, 'ctp[m]': -2, 'gtp[m]': -1, 'utp[m]': -5, 'ppi[m]': 9, 'trans2[m]': 1})
        self.assertEqual(len(self.model.observables), 4)
        self.assertEqual(self.model.observables.get_one(name='RNA Polymerase I observable in nucleus').id, 'obs_1')
        self.assertEqual(self.model.observables.get_one(name='RNA Polymerase I observable in nucleus').expression.expression,
            'complex1[n] + complex2[n]')
        self.assertEqual(self.model.observables.get_one(name='RNA Polymerase II observable in mitochondria').expression.expression,
            'complex3[m]')
        self.assertEqual(self.model.observables.get_one(name='RNA Polymerase III observable in mitochondria').expression.expression,
            'complex4[m]')
        self.assertEqual(self.model.observables.get_one(name='Combined RNA polymerase observable in mitochondria').expression.expression,
            'obs_2 + obs_3')

        # Test gen_rate_laws
        self.assertEqual(len(self.model.rate_laws), 3)
        self.assertEqual(self.model.rate_laws.get_one(id='transcription_trans1-forward').expression.expression,
            'k_cat_transcription_trans1 * obs_1 * '
            '(atp[n] / (atp[n] + K_m_transcription_trans1_atp * Avogadro * volume_n)) * '
            '(ctp[n] / (ctp[n] + K_m_transcription_trans1_ctp * Avogadro * volume_n)) * '
            '(gtp[n] / (gtp[n] + K_m_transcription_trans1_gtp * Avogadro * volume_n)) * '
            '(utp[n] / (utp[n] + K_m_transcription_trans1_utp * Avogadro * volume_n))')
        self.assertEqual(self.model.rate_laws.get_one(id='transcription_trans2-forward').expression.expression,
            'k_cat_transcription_trans2 * obs_2 * '
            '(atp[m] / (atp[m] + K_m_transcription_trans2_atp * Avogadro * volume_m)) * '
            '(ctp[m] / (ctp[m] + K_m_transcription_trans2_ctp * Avogadro * volume_m)) * '
            '(gtp[m] / (gtp[m] + K_m_transcription_trans2_gtp * Avogadro * volume_m)) * '
            '(utp[m] / (utp[m] + K_m_transcription_trans2_utp * Avogadro * volume_m))')
        self.assertEqual(self.model.rate_laws.get_one(id='transcription_trans3-forward').expression.expression,
            'k_cat_transcription_trans3 * obs_4 * '
            '(atp[m] / (atp[m] + K_m_transcription_trans3_atp * Avogadro * volume_m)) * '
            '(ctp[m] / (ctp[m] + K_m_transcription_trans3_ctp * Avogadro * volume_m)) * '
            '(gtp[m] / (gtp[m] + K_m_transcription_trans3_gtp * Avogadro * volume_m)) * '
            '(utp[m] / (utp[m] + K_m_transcription_trans3_utp * Avogadro * volume_m))')

        # Test calibrate_submodel
        self.assertEqual(self.model.parameters.get_one(id='K_m_transcription_trans1_utp').value, 1500/scipy.constants.Avogadro/5E-14)
        self.assertEqual(self.model.parameters.get_one(id='K_m_transcription_trans2_utp').value, 1500/scipy.constants.Avogadro/2.5E-14)
        self.assertEqual(self.model.parameters.get_one(id='k_cat_transcription_trans1').value, math.log(2)*(1/(20*3600) + 1/36000)*10/(0.5**4*200))
        self.assertEqual(self.model.parameters.get_one(id='k_cat_transcription_trans2').value, math.log(2)*(1/(20*3600) + 1/15000)*10/(0.5**4*100))
