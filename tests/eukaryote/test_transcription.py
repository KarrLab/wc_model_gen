""" Tests of transcription submodel generation

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-05-21
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import transcription
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

class TranscriptionSubmodelGeneratorTestCase(unittest.TestCase):

    def setUp(self):

        # Create KB content
        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>chr1\nATGCATGACTCTAGTTTAT\n'
                    '>chrM\nTTTatgaCTCTAGTTTACTTT\n')

        self.kb = wc_kb.KnowledgeBase()
        cell = self.kb.cell = wc_kb.Cell()

        nucleus = cell.compartments.create(id='n')
        mito = cell.compartments.create(id='m')

        chr1 = wc_kb.core.DnaSpeciesType(cell=cell, id='chr1', sequence_path=self.sequence_path)
        gene1 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene1', polymer=chr1, start=1, end=19)
        exon1 = wc_kb.eukaryote.GenericLocus(start=5, end=19)
        transcript1 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans1', 
            name='transcript1', gene=gene1, exons=[exon1])
        transcript1_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=transcript1, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        transcript1_spec = wc_kb.core.Species(species_type=transcript1, compartment=nucleus)
        transcript1_conc = wc_kb.core.Concentration(cell=cell, species=transcript1_spec, value=10.)

        chrM = wc_kb.core.DnaSpeciesType(cell=cell, id='chrM', sequence_path=self.sequence_path)
        gene2 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene2', polymer=chrM, start=1, end=18)
        exon2 = wc_kb.eukaryote.GenericLocus(start=1, end=10)
        transcript2 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans2', 
            name='transcript2', gene=gene2, exons=[exon2])
        transcript2_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=transcript2, 
            value='15000.0', value_type=wc_ontology['WC:float'])
        transcript2_spec = wc_kb.core.Species(species_type=transcript2, compartment=mito)
        transcript2_conc = wc_kb.core.Concentration(cell=cell, species=transcript2_spec, value=10.)

        gene3 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene3', polymer=chr1, start=1, end=19)
        exon3 = wc_kb.eukaryote.GenericLocus(start=1, end=15)
        transcript3 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans3', 
            name='transcript3', gene=gene3, exons=[exon3])
        transcript3_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=transcript3, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        transcript3_spec = wc_kb.core.Species(species_type=transcript3, compartment=nucleus)
        transcript3_conc = wc_kb.core.Concentration(cell=cell, species=transcript3_spec, value=10.)

        gene4 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene4', polymer=chr1, start=1, end=3)
        exon4 = wc_kb.eukaryote.GenericLocus(start=1, end=3)
        transcript4 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans4', 
            name='transcript4', gene=gene4, exons=[exon4])
        transcript4_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=transcript4, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        transcript4_spec = wc_kb.core.Species(species_type=transcript4, compartment=nucleus)
        transcript4_conc = wc_kb.core.Concentration(cell=cell, species=transcript4_spec, value=10.)

        gene5 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene5', polymer=chr1, start=1, end=3)
        exon5 = wc_kb.eukaryote.GenericLocus(start=1, end=3)
        transcript5 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans5', 
            name='transcript5', gene=gene5, exons=[exon5])
        transcript5_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=transcript5, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        transcript5_spec = wc_kb.core.Species(species_type=transcript5, compartment=nucleus)
        transcript5_conc = wc_kb.core.Concentration(cell=cell, species=transcript5_spec, value=0.)

        activator = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='activator')
        repressor = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='repressor')
        gene2_reg1 = gene2.regulatory_modules.create(
            transcription_factor_regulation=[wc_kb.eukaryote.TranscriptionFactorRegulation(
            transcription_factor=activator, direction=wc_kb.eukaryote.RegulatoryDirection.activation)])
        gene2_reg2 = gene2.regulatory_modules.create(
            transcription_factor_regulation=[wc_kb.eukaryote.TranscriptionFactorRegulation(
            transcription_factor=repressor, direction=wc_kb.eukaryote.RegulatoryDirection.repression)])                           

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

        for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType):
            model_species_type = model.species_types.create(id=i.id)
            model_compartment = model.compartments.get_one(id='m' if 'M' in i.gene.polymer.id else 'n')
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=model_species, 
                mean=10, units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()
        model.distribution_init_concentrations.get_one(id='dist-init-conc-trans5[n]').mean = 0.

        complexes = {'complex1': ('RNA Polymerase I','n'), 'complex2': ('RNA Polymerase II', 'n'), 
            'complex3': ('RNA Polymerase mitochondria', 'm'), 'complex4': ('RNA Polymerase III', 'n')}
        for k, v in complexes.items():
            model_species_type = model.species_types.create(id=k, name=v[0], 
                structure = wc_lang.ChemicalStructure(
                    empirical_formula = EmpiricalFormula('H'),
                    molecular_weight = 1.018,
                    charge = 1))
            model_compartment = model.compartments.get_one(id=v[1])
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=model_species, 
                mean=2500 if k=='complex2' else 100, units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()

        metabolic_participants = ['atp', 'ctp', 'gtp', 'utp', 'ppi', 'amp', 'cmp', 'gmp', 'ump', 'h2o', 'h', 'adp', 'pi']
        for i in metabolic_participants:
            model_species_type = model.species_types.create(id=i, name=i.upper())
            for c in ['n', 'm']:
                model_compartment = model.compartments.get_one(id=c)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
                conc_model = model.distribution_init_concentrations.create(species=model_species, 
                    mean=1500, units=unit_registry.parse_units('molecule'))
                conc_model.id = conc_model.gen_id()

        for i in ['activator', 'repressor']:
            model_species_type = model.species_types.create(id=i, name=i.upper())
            model_compartment = model.compartments.get_one(id='m')
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.create(species=model_species, 
                mean=3, units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()        

    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)
        gvar.transcript_ntp_usage = {}            

    def test_methods(self):           

        gen = transcription.TranscriptionSubmodelGenerator(self.kb, self.model, options={
            'rna_pol_pair': {'trans1': 'RNA Polymerase I', 'trans2': 'RNA Polymerase mitochondria', 
                            'trans3': 'RNA Polymerase II', 'trans4': 'RNA Polymerase II', 'trans5': 'RNA Polymerase III'},
            'polr_occupancy_width': 2,                
            })
        gen.run()

        model=self.model

        self.assertEqual(gvar.transcript_ntp_usage['trans1'], {'A': 4, 'U': 7, 'G': 2, 'C': 2, 'len': 15})

        # Test gen_reactions
        self.assertEqual([i.id for i in self.model.submodels], ['transcription'])

        # binding to non-specific site
        nuclear_genome_binding_site_conc = model.distribution_init_concentrations.get_one(
            id='dist-init-conc-polr_non_specific_binding_site[n]')
        mito_genome_binding_site_conc = model.distribution_init_concentrations.get_one(
            id='dist-init-conc-polr_non_specific_binding_site[m]')
        self.assertEqual(model.species_types.get_one(id='polr_non_specific_binding_site').structure.empirical_formula, EmpiricalFormula())
        self.assertEqual(model.species_types.get_one(id='polr_non_specific_binding_site').structure.molecular_weight, 0.)
        self.assertEqual(model.species_types.get_one(id='polr_non_specific_binding_site').structure.charge, 0)
        self.assertEqual(nuclear_genome_binding_site_conc.mean, 9)
        self.assertEqual(mito_genome_binding_site_conc.mean, 10)
        self.assertEqual(mito_genome_binding_site_conc.comments, 
            'Set to genome length divided by 2 bp to allow queueing of RNA polymerase during transcription')
        self.assertEqual(nuclear_genome_binding_site_conc.references[0].title, 
            'Structure and mechanism of the RNA Polymerase II transcription machinery')
        self.assertEqual(model.distribution_init_concentrations.get_one(
            id='dist-init-conc-complex1[n]').mean, 75)
        self.assertEqual(model.distribution_init_concentrations.get_one(
            id='dist-init-conc-complex1[n]').comments, 
            'The free pool is estimated to be three quarters of the total concentration')
        self.assertEqual(model.distribution_init_concentrations.get_one(
            id='dist-init-conc-complex1[n]').references[0].title, 
            'In vivo dynamics of RNA polymerase II transcription')
        self.assertEqual(model.species_types.get_one(id='complex1_bound_non_specific_site').structure.empirical_formula, EmpiricalFormula('H'))
        self.assertEqual(model.species_types.get_one(id='complex1_bound_non_specific_site').structure.molecular_weight, 1.018)
        self.assertEqual(model.species_types.get_one(id='complex1_bound_non_specific_site').structure.charge, 1)
        self.assertEqual(model.distribution_init_concentrations.get_one(
            id='dist-init-conc-complex1_bound_non_specific_site[n]').mean, 24)
        self.assertEqual(model.distribution_init_concentrations.get_one(
            id='dist-init-conc-complex3_bound_non_specific_site[m]').comments, 
            'Approximately 24.75 percent of RNA polymerase is bound to non-specific site')
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='non_specific_binding_complex1').participants},
            {'complex1[n]': -1, 'polr_non_specific_binding_site[n]': -1, 'complex1_bound_non_specific_site[n]': 1})
        self.assertEqual(model.reactions.get_one(id='non_specific_binding_complex2').submodel.id, 'transcription')

        # initiation
        self.assertEqual(model.species_types.get_one(id='gene1_binding_site').structure.empirical_formula, EmpiricalFormula())
        self.assertEqual(model.species_types.get_one(id='gene1_binding_site').structure.molecular_weight, 0.)
        self.assertEqual(model.species_types.get_one(id='gene1_binding_site').structure.charge, 0)
        self.assertEqual(model.species_types.get_one(id='complex1_bound_gene1').structure.empirical_formula, EmpiricalFormula('H'))
        self.assertEqual(model.species_types.get_one(id='complex1_bound_gene1').structure.molecular_weight, 1.018)
        self.assertEqual(model.species_types.get_one(id='complex1_bound_gene1').structure.charge, 1)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-gene1_binding_site[n]').mean, 10)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-gene2_binding_site[m]').mean, 10)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-gene4_binding_site[n]').mean, 2)
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='transcription_initiation_trans1').participants},
            {'gene1_binding_site[n]': -1, 'complex1_bound_non_specific_site[n]': -1, 'polr_non_specific_binding_site[n]': 1, 'complex1_bound_gene1[n]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='transcription_initiation_trans2').participants},
            {'gene2_binding_site[m]': -1, 'complex3_bound_non_specific_site[m]': -1, 'polr_non_specific_binding_site[m]': 1, 'complex3_bound_gene2[m]': 1})      
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='transcription_initiation_trans3').participants},
            {'gene3_binding_site[n]': -1, 'complex2_bound_non_specific_site[n]': -1, 'polr_non_specific_binding_site[n]': 1, 'complex2_bound_gene3[n]': 1,
            'atp[n]': -2, 'h2o[n]': -2, 'adp[n]': 2, 'pi[n]': 2, 'h[n]': 2})      
        
        # elongation        
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='transcription_elongation_trans1').participants}, 
            {'atp[n]': -5, 'ctp[n]': -3, 'gtp[n]': -3, 'utp[n]': -8, 'h2o[n]': -5,'ppi[n]': 19, 'trans1[n]': 1, 
            'amp[n]': 1, 'cmp[n]': 1, 'gmp[n]': 1, 'ump[n]': 1, 'h[n]': 5, 
            'complex1_bound_gene1[n]': -1, 'gene1_binding_site[n]': 1, 'complex1[n]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in model.reactions.get_one(id='transcription_elongation_trans2').participants}, 
            {'atp[m]': -4, 'ctp[m]': -3, 'gtp[m]': -2, 'utp[m]': -9, 'h2o[m]': -9, 'ppi[m]': 18, 'trans2[m]': 1, 
            'amp[m]': 2, 'cmp[m]': 1, 'gmp[m]': 1, 'ump[m]': 4,'h[m]': 9,
            'complex3_bound_gene2[m]': -1, 'gene2_binding_site[m]': 1, 'complex3[m]': 1})


        # Test gen_rate_laws
        self.assertEqual(len(model.rate_laws), 14)

        # binding to non-specific site
        self.assertEqual(model.rate_laws.get_one(id='non_specific_binding_complex1-forward').expression.expression,
            'k_non_specific_binding_complex1 * complex1[n]')
        self.assertEqual(model.rate_laws.get_one(id='non_specific_binding_complex2-forward').expression.expression,
            'k_non_specific_binding_complex2 * complex2[n]')
        self.assertEqual(model.rate_laws.get_one(id='non_specific_binding_complex3-forward').expression.expression,
            'k_non_specific_binding_complex3 * complex3[m]')

        # initiation
        self.assertEqual(model.observables.get_one(id='subtotal_complex1_n_1').expression.expression,
            'complex1_bound_gene1[n] + complex1[n] + complex1_bound_non_specific_site[n]')
        self.assertEqual(model.observables.get_one(id='total_complex1_n').expression.expression,
            'subtotal_complex1_n_1')
        self.assertEqual(model.parameters.get_one(id='total_nuclear_genome_binding').value, 9)
        self.assertEqual(model.parameters.get_one(id='total_nuclear_genome_binding').units, unit_registry.parse_units('molecule'))
        self.assertEqual(model.parameters.get_one(id='total_nuclear_genome_binding').comments, 'Set to genome length divided by 2 bp')
        self.assertEqual(model.parameters.get_one(id='total_nuclear_genome_binding').references[0].author, 'Steven Hahn')
        self.assertEqual(model.parameters.get_one(id='total_mitochondrial_genome_binding').value, 10)
        self.assertEqual(model.parameters.get_one(id='total_mitochondrial_genome_binding').units, unit_registry.parse_units('molecule'))
        self.assertEqual(model.parameters.get_one(id='total_mitochondrial_genome_binding').comments, 'Set to genome length divided by 2 bp')
        self.assertEqual(model.parameters.get_one(id='total_mitochondrial_genome_binding').references[0].pages, '394-403')
        self.assertEqual(model.functions.get_one(id='p_bound_1').expression.expression,
            '1 / (1 + total_nuclear_genome_binding / (total_complex1_n * 1) * exp(log(K_d_specific_polr / K_d_non_specific_polr)))')
        self.assertEqual(model.functions.get_one(id='p_bound_2').expression.expression,
            '1 / (1 + total_mitochondrial_genome_binding / (total_complex3_m * ((1 + activator[m] / (Ka_transcription_initiation_trans2_activator * '
            'Avogadro * volume_m) * f_transcription_initiation_trans2_activator) / (1 + activator[m] / '
            '(Ka_transcription_initiation_trans2_activator * Avogadro * volume_m))) * (1 / (1 + repressor[m] / '
            '(Kr_transcription_initiation_trans2_repressor * Avogadro * volume_m)))) * '
            'exp(log(K_d_specific_polr / K_d_non_specific_polr)))')
        self.assertEqual(model.rate_laws.get_one(id='transcription_initiation_trans1-forward').expression.expression,
            'p_bound_1 * k_specific_binding_complex1 * complex1_bound_non_specific_site[n] * max(min(gene1_binding_site[n] , 1) , 0)')
        self.assertEqual(model.rate_laws.get_one(id='transcription_initiation_trans2-forward').expression.expression,
            'p_bound_2 * k_specific_binding_complex3 * complex3_bound_non_specific_site[m] * max(min(gene2_binding_site[m] , 1) , 0)')

        # elongation
        self.assertEqual(model.rate_laws.get_one(id='transcription_elongation_trans2-forward').expression.expression,
            'k_cat_transcription_elongation_trans2 * complex3_bound_gene2[m] * '
            '(atp[m] / (atp[m] + K_m_transcription_elongation_trans2_atp * Avogadro * volume_m)) * '
            '(ctp[m] / (ctp[m] + K_m_transcription_elongation_trans2_ctp * Avogadro * volume_m)) * '
            '(gtp[m] / (gtp[m] + K_m_transcription_elongation_trans2_gtp * Avogadro * volume_m)) * '
            '(utp[m] / (utp[m] + K_m_transcription_elongation_trans2_utp * Avogadro * volume_m))')
        self.assertEqual(model.rate_laws.get_one(id='transcription_elongation_trans3-forward').expression.expression,
            'k_cat_transcription_elongation_trans3 * complex2_bound_gene3[n] * '
            '(atp[n] / (atp[n] + K_m_transcription_elongation_trans3_atp * Avogadro * volume_n)) * '
            '(ctp[n] / (ctp[n] + K_m_transcription_elongation_trans3_ctp * Avogadro * volume_n)) * '
            '(gtp[n] / (gtp[n] + K_m_transcription_elongation_trans3_gtp * Avogadro * volume_n)) * '
            '(utp[n] / (utp[n] + K_m_transcription_elongation_trans3_utp * Avogadro * volume_n))')

        
        # Test calibrate_submodel
        self.assertEqual(model.parameters.get_one(id='K_m_transcription_elongation_trans1_utp').value, 1500/scipy.constants.Avogadro/5E-14)
        self.assertEqual(model.parameters.get_one(id='K_m_transcription_elongation_trans1_utp').comments, 
            'The value was assumed to be 1.0 times the concentration of utp in nucleus')
        self.assertEqual(model.parameters.get_one(id='K_m_transcription_elongation_trans2_utp').value, 1500/scipy.constants.Avogadro/2.5E-14)
        self.assertEqual(model.parameters.get_one(id='Kr_transcription_initiation_trans2_repressor').value, 3/scipy.constants.Avogadro/2.5E-14)
        self.assertEqual(model.parameters.get_one(id='Kr_transcription_initiation_trans2_repressor').comments, 
            'The value was assumed to be 1.0 times the concentration of repressor in mitochondria')
        self.assertEqual(model.parameters.get_one(id='Ka_transcription_initiation_trans2_activator').value, 3/scipy.constants.Avogadro/2.5E-14)
        self.assertEqual(model.parameters.get_one(id='Ka_transcription_initiation_trans2_activator').comments, 
            'The value was assumed to be 1.0 times the concentration of activator in mitochondria')
        self.assertEqual(model.parameters.get_one(id='f_transcription_initiation_trans2_activator').value, 1.2)

        self.assertEqual(model.parameters.get_one(id='k_non_specific_binding_complex1').value, math.log(2)*(1/(20*3600) + 1/36000)*10/75)
        self.assertEqual(model.parameters.get_one(id='k_non_specific_binding_complex3').value, math.log(2)*(1/(20*3600) + 1/15000)*10/75)
        self.assertAlmostEqual(model.parameters.get_one(id='k_specific_binding_complex1').value, 
            math.log(2)*(1/(20*3600) + 1/36000)*10/(24*1/(1+9/100*math.exp(math.log(1e-09/1e-03)))), places=20)
        self.assertAlmostEqual(model.parameters.get_one(id='k_specific_binding_complex2').value, 
            2*math.log(2)*(1/(20*3600) + 1/36000)*10/(618*2/(1+9/2500*math.exp(math.log(1e-09/1e-03)))), places=20)
        self.assertAlmostEqual(model.parameters.get_one(id='k_specific_binding_complex3').value, 
            math.log(2)*(1/(20*3600) + 1/15000)*10/(24*1/(1+10/(100*((1+1.2)/(1+1))*(1/(1+1)))*math.exp(math.log(1e-09/1e-03)))), places=20)

        self.assertEqual(model.parameters.get_one(id='k_cat_transcription_elongation_trans1').value, math.log(2)*(1/(20*3600) + 1/36000)*10/(0.5**4*1))
        self.assertEqual(model.parameters.get_one(id='k_cat_transcription_elongation_trans2').value, math.log(2)*(1/(20*3600) + 1/15000)*10/(0.5**4*1))
        self.assertEqual(model.parameters.get_one(id='k_cat_transcription_elongation_trans4').value, math.log(2)*(1/(20*3600) + 1/36000)*10/(0.5**3*2))
        self.assertAlmostEqual(model.parameters.get_one(id='k_cat_transcription_elongation_trans5').value, 0.0029706307738283375, places=16)
        self.assertAlmostEqual(model.parameters.get_one(id='k_cat_transcription_elongation_trans5').comments, 
            'Set to the median value because it could not be determined from data')

    def test_global_vars(self):
        
        gvar.transcript_ntp_usage = {'trans1': {'A': 2, 'U': 2, 'G': 2, 'C': 2, 'len': 8}}
        gen = transcription.TranscriptionSubmodelGenerator(self.kb, self.model, options={
            'rna_pol_pair': {'trans1': 'RNA Polymerase I', 'trans2': 'RNA Polymerase mitochondria', 
                            'trans3': 'RNA Polymerase II', 'trans4': 'RNA Polymerase II', 'trans5': 'RNA Polymerase III'},
            'polr_occupancy_width': 2,                
            })
        gen.run()

        self.assertEqual(gvar.transcript_ntp_usage['trans2'], {'A': 2, 'U': 5, 'G': 1, 'C': 2, 'len': 10})

        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='transcription_elongation_trans1').participants}, 
            {'atp[n]': -5, 'ctp[n]': -3, 'gtp[n]': -3, 'utp[n]': -8, 'h2o[n]': -12,'ppi[n]': 19, 'trans1[n]': 1, 
            'amp[n]': 3, 'cmp[n]': 1, 'gmp[n]': 1, 'ump[n]': 6, 'h[n]': 12, 
            'complex1_bound_gene1[n]': -1, 'gene1_binding_site[n]': 1, 'complex1[n]': 1})
        self.assertEqual({i.species.id: i.coefficient for i in self.model.reactions.get_one(id='transcription_elongation_trans2').participants}, 
            {'atp[m]': -4, 'ctp[m]': -3, 'gtp[m]': -2, 'utp[m]': -9, 'h2o[m]': -9, 'ppi[m]': 18, 'trans2[m]': 1, 
            'amp[m]': 2, 'cmp[m]': 1, 'gmp[m]': 1, 'ump[m]': 4,'h[m]': 9,
            'complex3_bound_gene2[m]': -1, 'gene2_binding_site[m]': 1, 'complex3[m]': 1})
