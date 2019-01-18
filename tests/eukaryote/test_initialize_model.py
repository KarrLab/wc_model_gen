""" Tests of model initialization

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-11
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import core
from wc_model_gen.eukaryote import initialize_model
from wc_utils.util import chem
from wc_utils.util.ontology import wcm_ontology
import Bio.SeqUtils
import mendeleev
import os
import shutil
import tempfile
import unittest
import wc_kb
import wc_lang


class TestCase(unittest.TestCase):

    @staticmethod
    def set_options(options):
		
        option_dict = { 'gen_dna': False,
                        'gen_pre_rnas': False,
                        'gen_transcripts': False,
                        'gen_protein': False,
                        'gen_metabolites': False,
                        'gen_complexes': False,
                        'gen_distribution_init_concentrations': False,
                        'gen_observables': False,
                        'gen_kb_reactions': False,
                        'gen_kb_rate_laws': False,
                        }            
		
        for i in options:
            del option_dict[i] 

        return option_dict

    def setUp(self):

        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>chr1\nTTTATGAARGTNCTCATHAAYAARAAYGARCTCTAGTTTAT\n'
                    '>chrX\nATGCGTCA\n'
                    '>chrMT\nATGAARAARTTYCTCCTCACNCCNCTCTAATTT\n')
    
        self.kb = wc_kb.KnowledgeBase()
        cell = self.kb.cell = wc_kb.Cell()

        cell.properties.create(id='cell_volume', value=10400.)
        cell.properties.create(id='mean_doubling_time', value=20., units='hours')
        nucleus = cell.compartments.create(id='n', name='nucleus', volumetric_fraction=0.5)
        mito = cell.compartments.create(id='m', name='mitochondrion', volumetric_fraction=0.2)
        extra = cell.compartments.create(id='e', name='extracellular space')

        chr1 = wc_kb.core.DnaSpeciesType(cell=cell, id='chr1', name='chromosome 1', 
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene1 = wc_kb.eukaryote_schema.GeneLocus(polymer=chr1, start=1, end=36)
        rna1 = wc_kb.eukaryote_schema.PreRnaSpeciesType(cell=cell, id='rna1', name='rna1', gene=gene1)
        exon1 = wc_kb.eukaryote_schema.ExonLocus(polymer=chr1, start=4, end=36)
        transcript1 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans1', 
            name='transcript1', rna=rna1, exons=[exon1])
        transcript1_spec = wc_kb.core.Species(species_type=transcript1, compartment=nucleus)
        transcript1_conc = wc_kb.core.Concentration(cell=cell, species=transcript1_spec, value=0.02)
        cds1 = wc_kb.eukaryote_schema.CdsLocus(id='cds1', exon=exon1, start=4, end=36)        
        prot1 = wc_kb.eukaryote_schema.ProteinSpeciesType(cell=cell, id='prot1', name='protein1', 
            transcript=transcript1, coding_regions=[cds1])
        prot1_spec = wc_kb.core.Species(species_type=prot1, compartment=nucleus)
        prot1_conc = wc_kb.core.Concentration(cell=cell, species=prot1_spec, value=0.03)

        chrX = wc_kb.core.DnaSpeciesType(cell=cell, id='chrX', name='chromosome X', 
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene2 = wc_kb.eukaryote_schema.GeneLocus(polymer=chrX, start=1, end=4)
        rna2 = wc_kb.eukaryote_schema.PreRnaSpeciesType(cell=cell, id='rna2', name='rna2', gene=gene2)
        exon2 = wc_kb.eukaryote_schema.ExonLocus(polymer=chrX, start=1, end=4)
        transcript2 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans2',
            name='transcript2', rna=rna2, exons=[exon2])
        transcript2_spec = wc_kb.core.Species(species_type=transcript2, compartment=nucleus)
        transcript2_conc = wc_kb.core.Concentration(cell=cell, species=transcript2_spec, value=0.01)         

        chrMT = wc_kb.core.DnaSpeciesType(cell=cell, id='chrMT', name='mitochondrial chromosome', 
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene3 = wc_kb.eukaryote_schema.GeneLocus(polymer=chrMT, start=1, end=33)
        rna3 = wc_kb.eukaryote_schema.PreRnaSpeciesType(cell=cell, id='rna3', name='rna3', gene=gene3)
        exon3 = wc_kb.eukaryote_schema.ExonLocus(polymer=chrMT, start=1, end=30)
        transcript3 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans3', 
            name='transcript3',rna=rna3, exons=[exon3])
        transcript3_spec = wc_kb.core.Species(species_type=transcript3, compartment=mito)
        transcript3_conc = wc_kb.core.Concentration(cell=cell, species=transcript3_spec, value=0.05)
        cds3 = wc_kb.eukaryote_schema.CdsLocus(id='cds3', exon=exon3, start=1, end=30)        
        prot3 = wc_kb.eukaryote_schema.ProteinSpeciesType(cell=cell, id='prot3', name='protein3', 
            transcript=transcript3, coding_regions=[cds3])
        prot3_spec = wc_kb.core.Species(species_type=prot3, compartment=mito)
        prot3_conc = wc_kb.core.Concentration(cell=cell, species=prot3_spec, value=0.1)

        met1 = wc_kb.core.MetaboliteSpeciesType(cell=cell, id='met1', name='metabolite1', structure=(
            'InChI=1S'
            '/C10H14N5O7P'
            '/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)'
            '/p-2/t4-,6-,7-,10-'
            '/m1'
            '/s1'
        ))
        met1_spec1 = wc_kb.core.Species(species_type=met1, compartment=nucleus)
        met1_conc1 = wc_kb.core.Concentration(cell=cell, species=met1_spec1, value=0.1)
        met1_spec2 = wc_kb.core.Species(species_type=met1, compartment=extra)
        met1_conc2 = wc_kb.core.Concentration(cell=cell, species=met1_spec2, value=0.5)

        species_type_coeff1 = wc_kb.core.SpeciesTypeCoefficient(species_type=prot1, coefficient=2)
        species_type_coeff2 = wc_kb.core.SpeciesTypeCoefficient(species_type=prot3, coefficient=3)
        complex1 = wc_kb.core.ComplexSpeciesType(cell=cell, id='comp1', name='complex1',
            subunits=[species_type_coeff1, species_type_coeff2])

        expr1 = wc_kb.core.ObservableExpression(expression = '2.5 * prot1[n] + 1.3 * prot3[m]',
            species = [prot1_spec, prot3_spec])
        observable1 = wc_kb.core.Observable(cell=cell, id='obs1', name='observable1', expression=expr1)      

        expr2 = wc_kb.core.ObservableExpression(expression = '2.5 * prot1[n] / obs1',
            species = [prot1_spec], observables=[observable1])
        observable2 = wc_kb.core.Observable(cell=cell, id='obs2', name='observable2', expression=expr2)

    def tearDown(self):    
        shutil.rmtree(self.tmp_dirname)  

    def test_gen_compartments_and_parameters(self):        

        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([])}}).run()
        self.assertEqual(model.compartments.get_one(id='n').name, 'nucleus')
        self.assertEqual(model.compartments.get_one(id='n').mean_init_volume, 5200.)
        self.assertEqual(model.compartments.get_one(id='e').mean_init_volume, 10400E5)
        self.assertEqual(model.parameters.get_one(id='density_n').value, 1040.)
        self.assertEqual(model.parameters.get_one(id='density_e').value, 1000.)
        self.assertEqual(model.parameters.get_one(id='density_n').units, 'g l^-1')
        self.assertEqual(model.functions.get_one(id='volume_n').units, 'l')
        self.assertEqual(wc_lang.Function.expression.serialize(
            model.functions.get_one(id='volume_n').expression), 'n / density_n')
        self.assertEqual(wc_lang.Function.expression.serialize(
            model.functions.get_one(id='volume_e').expression), 'e / density_e')

        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').type, None)
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').units, 's')
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').value, 20*3600)

        self.kb.cell.properties.get_one(id='mean_doubling_time').units = 'minutes'
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([])}}).run()
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').units, 's')
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').value, 20*60)

        self.kb.cell.properties.get_one(id='mean_doubling_time').units = 'sec'
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([])}}).run()
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').units, 's')
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').value, 20)

    def test_gen_dna(self):        

        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_dna'])}}).run()

        chr1_model = model.species_types.get_one(id='chr1')
        chrX_model = model.species_types.get_one(id='chrX')
        chrMT_model = model.species_types.get_one(id='chrMT')

        self.assertEqual(chr1_model.name, 'chromosome 1')
        self.assertEqual(chr1_model.type, wcm_ontology['WCM:DNA'])
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=chr1_model)), True)
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=chrX_model)), True)
        self.assertEqual(all(i.compartment.id=='m' for i in model.species.get(species_type=chrMT_model)), True)
        
        dna = self.kb.cell.species_types.get_one(id='chrX')
        L = dna.get_len()
        self.assertEqual(chrX_model.empirical_formula, chem.EmpiricalFormula('C10H12N5O6P') * 2
                         + chem.EmpiricalFormula('C9H12N3O7P') * 2
                         + chem.EmpiricalFormula('C10H12N5O7P') * 2
                         + chem.EmpiricalFormula('C10H13N2O8P') * 2
                         - chem.EmpiricalFormula('OH') * (L - 1))

        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(dna.get_seq(),
                                            seq_type='DNA',
                                            circular=dna.circular,
                                            double_stranded=dna.double_stranded) \
            - 9 * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(chrX_model.molecular_weight, exp_mol_wt, places=0)
        self.assertEqual(chrX_model.charge, -L - 1)
        self.assertEqual(chrX_model.comments, '')

    def test_gen_pre_rnas(self):        
 
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_pre_rnas'])}}).run()

        rna1_model = model.species_types.get_one(id='rna1')
        rna2_model = model.species_types.get_one(id='rna2')
        rna3_model = model.species_types.get_one(id='rna3')

        self.assertEqual(rna1_model.name, 'rna1')
        self.assertEqual(rna2_model.type, wcm_ontology['WCM:pseudo_species'])
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=rna1_model)), True)
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=rna2_model)), True)
        self.assertEqual(all(i.compartment.id=='m' for i in model.species.get(species_type=rna3_model)), True)

        rna = self.kb.cell.species_types.get_one(id='rna2')
        L = rna.get_len()
        self.assertEqual(rna2_model.empirical_formula, chem.EmpiricalFormula('C10H12N5O7P') * 1
                         + chem.EmpiricalFormula('C9H12N3O8P') * 1
                         + chem.EmpiricalFormula('C10H12N5O8P') * 1
                         + chem.EmpiricalFormula('C9H11N2O9P') * 1
                         - chem.EmpiricalFormula('OH') * (L - 1))

        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(rna.get_seq()) \
            - (L + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(rna2_model.molecular_weight, exp_mol_wt, places=0)
        self.assertEqual(rna2_model.charge, -L - 1)
        self.assertEqual(rna2_model.comments, '')
        
    def test_gen_transcripts(self):
    
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_transcripts'])}}).run()

        transcript1_model = model.species_types.get_one(id='trans1')
        transcript2_model = model.species_types.get_one(id='trans2')
        transcript3_model = model.species_types.get_one(id='trans3')

        self.assertEqual(transcript1_model.name, 'transcript1')
        self.assertEqual(transcript3_model.type, wcm_ontology['WCM:RNA'])
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=transcript1_model)), True)
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=transcript2_model)), True)
        self.assertEqual(all(i.compartment.id=='m' for i in model.species.get(species_type=transcript3_model)), True)

        transcript = self.kb.cell.species_types.get_one(id='trans2')
        L = transcript.get_len()
        self.assertEqual(transcript2_model.empirical_formula, chem.EmpiricalFormula('C10H12N5O7P') * 1
                         + chem.EmpiricalFormula('C9H12N3O8P') * 1
                         + chem.EmpiricalFormula('C10H12N5O8P') * 1
                         + chem.EmpiricalFormula('C9H11N2O9P') * 1
                         - chem.EmpiricalFormula('OH') * (L - 1))

        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript.get_seq()) \
            - (L + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(transcript2_model.molecular_weight, exp_mol_wt, places=0)
        self.assertEqual(transcript2_model.charge, -L - 1)
        self.assertEqual(transcript2_model.comments, '')

    def test_gen_protein(self):

        self.kb.cell.species_types.get_one(id='prot3').comments = 'Testing'
    
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_protein'])}}).run()

        prot1_model = model.species_types.get_one(id='prot1')
        prot3_model = model.species_types.get_one(id='prot3')

        self.assertEqual(prot1_model.name, 'protein1')
        self.assertEqual(prot3_model.type, wcm_ontology['WCM:protein'])
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=prot1_model)), True)
        self.assertEqual(all(i.compartment.id=='m' for i in model.species.get(species_type=prot3_model)), True)
        self.assertEqual(prot1_model.empirical_formula, chem.EmpiricalFormula('C53H96N14O15S1'))
        self.assertAlmostEqual(prot1_model.molecular_weight, 1201.49, delta=0.3)
        self.assertEqual(prot1_model.charge, 1)
        self.assertEqual(prot1_model.comments, '')
        self.assertEqual(prot3_model.comments, 'Testing')

    def test_gen_metabolites(self):
    
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_metabolites'])}}).run()

        met1_model = model.species_types.get_one(id='met1')
        
        self.assertEqual(met1_model.name, 'metabolite1')
        self.assertEqual(met1_model.type, wcm_ontology['WCM:metabolite'])
        self.assertEqual(met1_model.structure, 'InChI=1S'
            '/C10H14N5O7P'
            '/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)'
            '/p-2/t4-,6-,7-,10-'
            '/m1'
            '/s1')
        self.assertEqual(set([i.compartment.id for i in model.species.get(species_type=met1_model)]), set(['n', 'e']))
        self.assertEqual(met1_model.empirical_formula, chem.EmpiricalFormula('C10H12N5O7P'))
        self.assertAlmostEqual(met1_model.molecular_weight, 345.20530, places=4)
        self.assertEqual(met1_model.charge, -2)
        self.assertEqual(met1_model.comments, '')

    def test_complexes(self):        
    
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_complexes'])}}).run()

        comp1_model = model.species_types.get_one(id='comp1')
        
        self.assertEqual(comp1_model.name, 'complex1')
        self.assertEqual(comp1_model.type, wcm_ontology['WCM:pseudo_species'])
        self.assertEqual(set([i.compartment.id for i in model.species.get(species_type=comp1_model)]), set(['n', 'm']))
        self.assertEqual(comp1_model.empirical_formula, chem.EmpiricalFormula('C53H96N14O15S1') * 2
                            + chem.EmpiricalFormula('C53H91N11O11S1') * 3)
        self.assertAlmostEqual(comp1_model.molecular_weight, 5674.27, delta=0.3)
        self.assertEqual(comp1_model.charge, 8)
        self.assertEqual(comp1_model.comments, '')

    def test_gen_distribution_init_concentrations(self):

        test_conc = self.kb.cell.concentrations.get_one(value=0.5)
        test_conc.comments = 'Testing'
        test_conc.references.append(wc_kb.core.Reference(id='ref1', standard_id='doi:1234'))
        test_conc.database_references.append(wc_kb.core.DatabaseReference(database='ECMDB', id='12345'))

        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_distribution_init_concentrations'])}}).run()
        
        met1_nucleus = model.distribution_init_concentrations.get_one(id='dist-init-conc-met1[n]')
        met1_extra = model.distribution_init_concentrations.get_one(id='dist-init-conc-met1[e]')
        
        self.assertEqual(met1_nucleus.species, model.species.get_one(id='met1[n]'))
        self.assertEqual(met1_extra.mean, 0.5)
        self.assertEqual(met1_nucleus.units, wc_lang.ConcentrationUnit.M)
        self.assertEqual(met1_nucleus.comments, '')
        self.assertEqual(met1_extra.comments, 'Testing')
        self.assertEqual(met1_nucleus.references, [])
        self.assertEqual(met1_extra.references[0].id, 'ref1')
        self.assertEqual(met1_extra.references[0].name, 'doi:1234')
        self.assertEqual(met1_nucleus.db_refs, [])
        self.assertEqual(met1_extra.db_refs[0].serialize(), 'ECMDB: 12345')

    def test_gen_observables(self):

        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel], 
            options={'component': {'InitializeModel': self.set_options(['gen_protein', 'gen_observables'])}}).run()

        self.assertEqual(model.observables.get_one(id='obs1').name, 'observable1')
        self.assertEqual(model.observables.get_one(id='obs1').expression.expression, '2.5 * prot1[n] + 1.3 * prot3[m]')
        self.assertEqual(set([i.species_type.id for i in model.observables.get_one(id='obs1').expression.species]), 
            set(['prot1', 'prot3']))
        self.assertEqual(set([i.compartment.id for i in model.observables.get_one(id='obs1').expression.species]), 
            set(['n', 'm']))    
        self.assertEqual(model.observables.get_one(id='obs2').expression.expression, '2.5 * prot1[n] / obs1')
        self.assertEqual(set([i.species_type.id for i in model.observables.get_one(id='obs2').expression.species]), 
            set(['prot1']))
        self.assertEqual(set([i.compartment.id for i in model.observables.get_one(id='obs2').expression.species]), 
            set(['n']))
        self.assertEqual(model.observables.get_one(id='obs2').expression.observables[0],  
            model.observables.get_one(id='obs1'))     

    def test_gen_kb_reactions(self):
        pass 

    def test_gen_kb_rate_laws(self):
        pass         

        