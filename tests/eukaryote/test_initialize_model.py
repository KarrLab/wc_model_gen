""" Tests of model initialization

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-11
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import core
from wc_model_gen.eukaryote import initialize_model
from wc_utils.util import chem
from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
from wc_lang.core import ChemicalStructure
import Bio.SeqUtils
import mendeleev
import os
import scipy.constants
import shutil
import tempfile
import unittest
import wc_kb
import wc_lang


class TestCase(unittest.TestCase):

    @staticmethod
    def set_options(options):

        option_dict = { 'gen_dna': False,
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
                    '>chrM\nATGAARAARTTYCTCCTCACNCCNCTCTAATTT\n')
    
        self.kb = wc_kb.KnowledgeBase(id='test_kb', version='0.0.1')
        cell = self.kb.cell = wc_kb.Cell(id='test_cell')

        cell.taxon = 9606

        ref = wc_kb.core.Reference(cell=cell, id='ref', authors='John Smith', year=2018, comments='No comment')
        cell.parameters.create(id='cell_volume', value=10400., references=[ref])
        cell.parameters.create(id='mean_doubling_time', value=20., units=unit_registry.parse_units('hour'))        

        nucleus = cell.compartments.create(id='n', name='nucleus', volumetric_fraction=0.5)
        nucleus_membrane = cell.compartments.create(id='n_m', name='nucleus membrane')
        mito = cell.compartments.create(id='m', name='mitochondrion', volumetric_fraction=0.2)
        extra = cell.compartments.create(id='e', name='extracellular space')

        chr1 = wc_kb.core.DnaSpeciesType(cell=cell, id='chr1', name='chromosome 1', ploidy=2,
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene1 = wc_kb.eukaryote_schema.GeneLocus(cell=cell, id='gene1', polymer=chr1, start=1, end=36)
        exon1 = wc_kb.eukaryote_schema.GenericLocus(start=4, end=36)
        transcript1 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans1',
            name='transcript1', gene=gene1, exons=[exon1])
        transcript1_spec = wc_kb.core.Species(species_type=transcript1, compartment=nucleus)
        transcript1_conc = wc_kb.core.Concentration(cell=cell, species=transcript1_spec, value=0.02)
        transcript1_conc.id = transcript1_conc.serialize()
        cds1 = wc_kb.eukaryote_schema.GenericLocus(start=4, end=36)        
        prot1 = wc_kb.eukaryote_schema.ProteinSpeciesType(cell=cell, id='prot1', name='protein1', 
            transcript=transcript1, coding_regions=[cds1])
        prot1_spec = wc_kb.core.Species(species_type=prot1, compartment=nucleus)
        prot1_conc = wc_kb.core.Concentration(cell=cell, species=prot1_spec, value=0.03)
        prot1_conc.id = prot1_conc.serialize()

        chrX = wc_kb.core.DnaSpeciesType(cell=cell, id='chrX', name='chromosome X', ploidy=1, 
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene2 = wc_kb.eukaryote_schema.GeneLocus(cell=cell, id='gene2', polymer=chrX, start=1, end=4)
        exon2 = wc_kb.eukaryote_schema.GenericLocus(start=1, end=4)
        transcript2 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans2',
            name='transcript2', gene=gene2, exons=[exon2])
        transcript2_spec = wc_kb.core.Species(species_type=transcript2, compartment=nucleus)
        transcript2_conc = wc_kb.core.Concentration(cell=cell, species=transcript2_spec, value=0.01)
        transcript2_conc.id = transcript2_conc.serialize()         

        chrM = wc_kb.core.DnaSpeciesType(cell=cell, id='chrM', name='mitochondrial chromosome', ploidy=150,
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene3 = wc_kb.eukaryote_schema.GeneLocus(cell=cell, id='gene3', polymer=chrM, start=1, end=33)
        exon3 = wc_kb.eukaryote_schema.GenericLocus(start=1, end=30)
        transcript3 = wc_kb.eukaryote_schema.TranscriptSpeciesType(cell=cell, id='trans3',
            name='transcript3', gene=gene3, exons=[exon3])
        transcript3_spec = wc_kb.core.Species(species_type=transcript3, compartment=mito)
        transcript3_conc = wc_kb.core.Concentration(cell=cell, species=transcript3_spec, value=0.05)
        transcript3_conc.id = transcript3_conc.serialize()
        cds3 = wc_kb.eukaryote_schema.GenericLocus(start=1, end=30)        
        prot3 = wc_kb.eukaryote_schema.ProteinSpeciesType(cell=cell, id='prot3', name='protein3', 
            transcript=transcript3, coding_regions=[cds3])
        prot3_spec = wc_kb.core.Species(species_type=prot3, compartment=mito)
        prot3_conc = wc_kb.core.Concentration(cell=cell, species=prot3_spec, value=0.1)
        prot3_conc.id = prot3_conc.serialize()

        met1 = wc_kb.core.MetaboliteSpeciesType(cell=cell, id='met1', name='metabolite1')
        met1_structure = wc_kb.core.SpeciesTypeProperty(property='structure', species_type=met1, value=
            'InChI=1S'
            '/C10H14N5O7P'
            '/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)'
            '/p-2/t4-,6-,7-,10-'
            '/m1'
            '/s1',
            value_type=wc_ontology['WC:string']
            )
        met1_structure.id = met1_structure.gen_id()

        met1_spec1 = wc_kb.core.Species(species_type=met1, compartment=nucleus)
        met1_conc1 = wc_kb.core.Concentration(cell=cell, species=met1_spec1, value=0.5)
        met1_conc1.id = met1_conc1.serialize()
        met1_spec2 = wc_kb.core.Species(species_type=met1, compartment=extra)

        species_type_coeff1 = wc_kb.core.SpeciesTypeCoefficient(species_type=prot1, coefficient=2)
        species_type_coeff2 = wc_kb.core.SpeciesTypeCoefficient(species_type=met1, coefficient=3)
        complex1 = wc_kb.core.ComplexSpeciesType(cell=cell, id='comp1', name='complex1',
            subunits=[species_type_coeff1, species_type_coeff2])

        species_type_coeff3 = wc_kb.core.SpeciesTypeCoefficient(species_type=prot1, coefficient=6)
        species_type_coeff4 = wc_kb.core.SpeciesTypeCoefficient(species_type=prot3, coefficient=4)
        complex2 = wc_kb.core.ComplexSpeciesType(cell=cell, id='comp2', name='complex2',
            subunits=[species_type_coeff3, species_type_coeff4])
        complex2_spec1 = wc_kb.core.Species(species_type=complex2, compartment=nucleus)
        complex2_conc1 = wc_kb.core.Concentration(cell=cell, species=complex2_spec1, value=0.)
        complex2_conc1.id = complex2_conc1.serialize()
        complex2_spec2 = wc_kb.core.Species(species_type=complex2, compartment=extra)
        complex2_conc2 = wc_kb.core.Concentration(cell=cell, species=complex2_spec2, value=0.)
        complex2_conc2.id = complex2_conc2.serialize()

        expr1 = wc_kb.core.ObservableExpression(expression = '2.5 * prot1[n] + 1.3 * prot3[m]',
            species = [prot1_spec, prot3_spec])
        observable1 = wc_kb.core.Observable(cell=cell, id='obs1', name='observable1', expression=expr1)

        expr2 = wc_kb.core.ObservableExpression(expression = '2.5 * prot1[n] / obs1',
            species = [prot1_spec], observables=[observable1])
        observable2 = wc_kb.core.Observable(cell=cell, id='obs2', name='observable2', expression=expr2)

        participant1 = wc_kb.core.SpeciesCoefficient(species=met1_spec1, coefficient=1)
        participant2 = wc_kb.core.SpeciesCoefficient(species=met1_spec2, coefficient=-1)
        reaction1 = wc_kb.core.Reaction(cell=cell, id='r1', name='reaction1',
            participants=[participant1, participant2], reversible=True)

        identifier = wc_kb.core.Identifier(namespace='Sabio-RK', id='1234')
        kcat_r1_forward = wc_kb.core.Parameter(cell=cell, id='k_cat_r1_forward', value=0.2, 
            units=unit_registry.parse_units('s^-1'), identifiers=[identifier])
        kcat_r1_backward = wc_kb.core.Parameter(cell=cell, id='k_cat_r1_backward', value=0.02, 
            units=unit_registry.parse_units('s^-1'))
        Km_r1_backward = wc_kb.core.Parameter(cell=cell, id='K_m_r1_backward_met1_n', value=0.3,
            units=unit_registry.parse_units('M'))
        forward_exp = wc_kb.core.RateLawExpression(
            expression='k_cat_r1_forward * comp2[e] * met1[e]',
            parameters=[kcat_r1_forward],
            species=[complex2_spec2, met1_spec2])
        backward_exp = wc_kb.core.RateLawExpression(
            expression='k_cat_r1_backward * comp2[n] * met1[n] * obs2 / (K_m_r1_backward_met1_n + met1[n])',
            parameters=[kcat_r1_backward, Km_r1_backward],
            species=[complex2_spec1, met1_spec1],
            observables=[observable2])
        forward_rate_law = wc_kb.core.RateLaw(
            reaction=reaction1,
            direction=wc_kb.core.RateLawDirection.forward,
            expression=forward_exp)
        forward_rate_law.id = forward_rate_law.gen_id()
        backward_rate_law = wc_kb.core.RateLaw(
            reaction=reaction1,
            direction=wc_kb.core.RateLawDirection.backward,
            expression=backward_exp)
        backward_rate_law.id = backward_rate_law.gen_id()            
        
    def tearDown(self):    
        shutil.rmtree(self.tmp_dirname)  

    def test_gen_taxon_compartments_parameters(self):

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([])}}).run()

        self.assertEqual(model.taxon.id, '9606')
        self.assertEqual(model.taxon.name, 'Homo sapiens')
        self.assertEqual(model.taxon.rank, wc_lang.TaxonRank.species)

        self.assertEqual(model.parameters.get_one(id='Avogadro').value, scipy.constants.Avogadro)
        self.assertEqual(model.parameters.get_one(id='Avogadro').type, None)
        self.assertEqual(model.parameters.get_one(id='Avogadro').units, unit_registry.parse_units('molecule mol^-1'))

        self.assertEqual(model.compartments.get_one(id='n').name, 'nucleus')
        self.assertAlmostEqual(model.compartments.get_one(id='n').init_volume.mean, 5199.999998548484)
        self.assertAlmostEqual(model.compartments.get_one(id='n_m').init_volume.mean, 1.4515160909356243e-06)
        self.assertEqual(model.compartments.get_one(id='e').init_volume.mean, 1.0)

        self.assertEqual(model.parameters.get_one(id='density_n').value, 1040.)
        self.assertEqual(model.parameters.get_one(id='density_n_m').value, 1160.)
        self.assertEqual(model.parameters.get_one(id='density_e').value, 1000.)
        self.assertEqual(model.parameters.get_one(id='density_n').units, unit_registry.parse_units('g l^-1'))

        self.assertEqual(model.functions.get_one(id='volume_n').units, unit_registry.parse_units('l'))
        self.assertEqual(wc_lang.Function.expression.serialize(
            model.functions.get_one(id='volume_n').expression), 'n / density_n')
        self.assertEqual(wc_lang.Function.expression.serialize(
            model.functions.get_one(id='volume_n_m').expression), 'n_m / density_n_m')
        
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').value, 0.2)
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').units, unit_registry.parse_units('s^-1'))

        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').type, wc_ontology['WC:k_cat'])
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').identifiers[0].namespace, 'Sabio-RK')
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').identifiers[0].id, '1234')
        self.assertEqual(model.parameters.get_one(id='K_m_r1_backward_met1_n').value, 0.3)
        self.assertEqual(model.parameters.get_one(id='K_m_r1_backward_met1_n').units, unit_registry.parse_units('M'))
        self.assertEqual(model.parameters.get_one(id='K_m_r1_backward_met1_n').type, wc_ontology['WC:K_m'])

        self.assertEqual(model.parameters.get_one(id='cell_volume').references[0].id, 'ref')
        self.assertEqual(model.parameters.get_one(id='cell_volume').references[0].author, 'John Smith')
        self.assertEqual(model.parameters.get_one(id='cell_volume').references[0].year, 2018)
        self.assertEqual(model.parameters.get_one(id='cell_volume').references[0].comments, 'No comment')

        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').type, None)
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').units, unit_registry.parse_units('s'))
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').value, 20*3600)

        self.kb.cell.parameters.get_one(id='mean_doubling_time').units = unit_registry.parse_units('minute')
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([])}}).run()
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').units, unit_registry.parse_units('s'))
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').value, 20*60)

        self.kb.cell.parameters.get_one(id='mean_doubling_time').units = unit_registry.parse_units('s')
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([])}}).run()
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').units, unit_registry.parse_units('s'))
        self.assertEqual(model.parameters.get_one(id='mean_doubling_time').value, 20)

    def test_gen_dna(self):

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_dna'])}}).run()

        chr1_model = model.species_types.get_one(id='chr1')
        chrX_model = model.species_types.get_one(id='chrX')
        chrM_model = model.species_types.get_one(id='chrM')

        self.assertEqual(chr1_model.name, 'chromosome 1')
        self.assertEqual(chr1_model.type, wc_ontology['WC:DNA'])
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=chr1_model)), True)
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=chrX_model)), True)
        self.assertEqual(all(i.compartment.id=='m' for i in model.species.get(species_type=chrM_model)), True)
        
        dna = self.kb.cell.species_types.get_one(id='chrX')
        L = dna.get_len()
        self.assertEqual(chrX_model.structure.empirical_formula, chem.EmpiricalFormula('C10H12N5O6P') * 2
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
        self.assertAlmostEqual(chrX_model.structure.molecular_weight, exp_mol_wt, places=0)
        self.assertEqual(chrX_model.structure.charge, -L - 1)
        self.assertEqual(chrX_model.comments, '')

    def test_gen_transcripts(self):

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_transcripts'])}}).run()

        transcript1_model = model.species_types.get_one(id='trans1')
        transcript2_model = model.species_types.get_one(id='trans2')
        transcript3_model = model.species_types.get_one(id='trans3')

        self.assertEqual(transcript1_model.name, 'transcript1')
        self.assertEqual(transcript3_model.type, wc_ontology['WC:RNA'])
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=transcript1_model)), True)
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=transcript2_model)), True)
        self.assertEqual(all(i.compartment.id=='m' for i in model.species.get(species_type=transcript3_model)), True)

        transcript = self.kb.cell.species_types.get_one(id='trans2')
        L = transcript.get_len()
        self.assertEqual(transcript2_model.structure.empirical_formula, chem.EmpiricalFormula('C10H12N5O7P') * 1
                         + chem.EmpiricalFormula('C9H12N3O8P') * 1
                         + chem.EmpiricalFormula('C10H12N5O8P') * 1
                         + chem.EmpiricalFormula('C9H11N2O9P') * 1
                         - chem.EmpiricalFormula('OH') * (L - 1))

        exp_mol_wt = \
            + Bio.SeqUtils.molecular_weight(transcript.get_seq()) \
            - (L + 1) * mendeleev.element('H').atomic_weight
        self.assertAlmostEqual(transcript2_model.structure.molecular_weight, exp_mol_wt, places=0)
        self.assertEqual(transcript2_model.structure.charge, -L - 1)
        self.assertEqual(transcript2_model.comments, '')

    def test_gen_protein(self):

        self.kb.cell.species_types.get_one(id='prot3').comments = 'Testing'

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_protein'])}}).run()

        prot1_model = model.species_types.get_one(id='prot1')
        prot3_model = model.species_types.get_one(id='prot3')

        self.assertEqual(prot1_model.name, 'protein1')
        self.assertEqual(prot3_model.type, wc_ontology['WC:protein'])
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=prot1_model)), True)
        self.assertEqual(all(i.compartment.id=='m' for i in model.species.get(species_type=prot3_model)), True)
        self.assertEqual(prot1_model.structure.empirical_formula, chem.EmpiricalFormula('C53H96N14O15S1'))
        self.assertAlmostEqual(prot1_model.structure.molecular_weight, 1201.49, delta=0.3)
        self.assertEqual(prot1_model.structure.charge, 1)
        self.assertEqual(prot1_model.comments, '')
        self.assertEqual(prot3_model.comments, 'Testing')

    def test_gen_metabolites(self):

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_metabolites'])}}).run()

        met1_model = model.species_types.get_one(id='met1')

        self.assertEqual(met1_model.name, 'metabolite1')
        self.assertEqual(met1_model.type, wc_ontology['WC:metabolite'])
        self.assertEqual(met1_model.structure.value, 'InChI=1S'
            '/C10H14N5O7P'
            '/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)'
            '/p-2/t4-,6-,7-,10-'
            '/m1'
            '/s1')
        self.assertEqual(set([i.compartment.id for i in model.species.get(species_type=met1_model)]), set(['n', 'e']))
        self.assertEqual(met1_model.structure.empirical_formula, chem.EmpiricalFormula('C10H12N5O7P'))
        self.assertAlmostEqual(met1_model.structure.molecular_weight, 345.20530, places=4)
        self.assertEqual(met1_model.structure.charge, -2)
        self.assertEqual(met1_model.comments, '')

    def test_complexes(self):

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_complexes'])}}).run()

        comp1_model = model.species_types.get_one(id='comp1')
        self.assertEqual(comp1_model.name, 'complex1')
        self.assertEqual(comp1_model.type, wc_ontology['WC:pseudo_species'])
        self.assertEqual(set([i.compartment.id for i in model.species.get(species_type=comp1_model)]), set(['n']))
        self.assertEqual(comp1_model.structure.empirical_formula, chem.EmpiricalFormula('C53H96N14O15S1') * 2
                            + chem.EmpiricalFormula('C10H12N5O7P') * 3)
        self.assertAlmostEqual(comp1_model.structure.molecular_weight, 3438.5959, delta=0.3)
        self.assertEqual(comp1_model.structure.charge, -4)
        self.assertEqual(comp1_model.comments, '')

        comp2_model = model.species_types.get_one(id='comp2')
        self.assertEqual(set([i.compartment.id for i in model.species.get(species_type=comp2_model)]), set(['n', 'e']))

    def test_gen_distribution_init_concentrations(self):

        test_conc = self.kb.cell.concentrations.get_one(value=0.5)
        test_conc.comments = 'Testing'
        test_conc.references.append(wc_kb.core.Reference(id='ref1', title='Title1', authors='Author1', 
            journal='Journal1', volume='1', issue='1', pages='20', year=1999, comments='xyz'))
        test_conc.identifiers.append(wc_kb.core.Identifier(namespace='ECMDB', id='12345'))

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([
                'gen_dna', 'gen_distribution_init_concentrations'])}}).run()

        met1_nucleus = model.distribution_init_concentrations.get_one(id='dist-init-conc-met1[n]')

        self.assertEqual(met1_nucleus.species, model.species.get_one(id='met1[n]'))
        self.assertEqual(met1_nucleus.mean, 0.5*scipy.constants.Avogadro*(0.5*10400.-4.836E-09*(0.5*10400.)**(2/3)))
        self.assertEqual(met1_nucleus.units, unit_registry.parse_units('molecule'))
        self.assertEqual(met1_nucleus.comments, 'Testing')
        self.assertEqual(met1_nucleus.references[0].id, 'ref1')
        self.assertEqual(met1_nucleus.references[0].title, 'Title1')
        self.assertEqual(met1_nucleus.references[0].author, 'Author1')
        self.assertEqual(met1_nucleus.references[0].publication, 'Journal1')
        self.assertEqual(met1_nucleus.references[0].volume, '1')
        self.assertEqual(met1_nucleus.references[0].issue, '1')
        self.assertEqual(met1_nucleus.references[0].pages, '20')
        self.assertEqual(met1_nucleus.references[0].year, 1999)
        self.assertEqual(met1_nucleus.references[0].comments, 'xyz')
        self.assertEqual(met1_nucleus.references[0].type, wc_ontology['WC:article'])
        self.assertEqual(met1_nucleus.identifiers[0].serialize(), 'ECMDB: 12345')

        test_conc.units = unit_registry.parse_units('molecule')

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([
                'gen_dna', 'gen_distribution_init_concentrations'])}}).run()
        
        met1_nucleus = model.distribution_init_concentrations.get_one(id='dist-init-conc-met1[n]')
        self.assertEqual(met1_nucleus.mean, 0.5)
        self.assertEqual(met1_nucleus.units, unit_registry.parse_units('molecule'))

        chr1_conc = model.distribution_init_concentrations.get_one(id='dist-init-conc-chr1[n]')
        chrX_conc = model.distribution_init_concentrations.get_one(id='dist-init-conc-chrX[n]')
        chrM_conc = model.distribution_init_concentrations.get_one(id='dist-init-conc-chrM[m]')

        self.assertEqual(chr1_conc.mean, 2)
        self.assertEqual(chrX_conc.mean, 1)
        self.assertEqual(chrM_conc.mean, 150)

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

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([
                'gen_protein', 'gen_metabolites', 'gen_complexes', 'gen_kb_reactions'])}}).run()

        self.assertEqual(len(model.reactions), 1)
        self.assertEqual(len(model.submodels), 1)

        self.assertEqual(model.reactions.get_one(id='r1_kb').name, 'reaction1')
        self.assertEqual(model.reactions.get_one(id='r1_kb').submodel.id, 'Metabolism')
        self.assertEqual(model.reactions.get_one(id='r1_kb').reversible, True)
        self.assertEqual(model.reactions.get_one(id='r1_kb').comments, '')
        self.assertEqual([(i.species.id, i.coefficient) for i in model.reactions.get_one(id='r1_kb').participants],
            [('met1[n]', 1), ('met1[e]', -1)])
        
    def test_gen_kb_rate_laws(self):
        
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel], 
            options={'component': {'InitializeModel': self.set_options(['gen_dna', 
                'gen_protein', 'gen_metabolites', 'gen_complexes', 'gen_observables', 
                'gen_distribution_init_concentrations', 'gen_kb_reactions', 
                'gen_kb_rate_laws'])}}).run()

        self.assertEqual(len(model.rate_laws), 2)

        self.assertEqual(model.rate_laws.get_one(id='r1_kb-forward').expression.expression,
            'k_cat_r1_forward * comp2[e] * met1[e]')
        self.assertEqual(model.rate_laws.get_one(id='r1_kb-forward').reaction, model.reactions.get_one(id='r1_kb'))
        self.assertEqual(model.rate_laws.get_one(id='r1_kb-forward').direction, wc_lang.RateLawDirection.forward)
        self.assertEqual(model.rate_laws.get_one(id='r1_kb-forward').comments, '')
        self.assertEqual([i.id for i in model.rate_laws.get_one(id='r1_kb-forward').expression.parameters],
            ['k_cat_r1_forward'])
        self.assertEqual(set([i.id for i in model.rate_laws.get_one(id='r1_kb-forward').expression.species]),
            set(['comp2[e]', 'met1[e]']))
        self.assertEqual(model.rate_laws.get_one(id='r1_kb-forward').expression.observables, [])

        self.assertEqual(model.rate_laws.get_one(id='r1_kb-backward').expression.expression,
            'k_cat_r1_backward * comp2[n] * met1[n] * obs2 / '
            '(K_m_r1_backward_met1_n * Avogadro * volume_n + met1[n])')
        self.assertEqual(model.rate_laws.get_one(id='r1_kb-backward').reaction, model.reactions.get_one(id='r1_kb'))
        self.assertEqual(model.rate_laws.get_one(id='r1_kb-backward').direction, wc_lang.RateLawDirection.backward)
        self.assertEqual(set([i.id for i in model.rate_laws.get_one(id='r1_kb-backward').expression.parameters]),
            set(['k_cat_r1_backward', 'K_m_r1_backward_met1_n', 'Avogadro']))
        self.assertEqual(set([i.id for i in model.rate_laws.get_one(id='r1_kb-backward').expression.species]),
            set(['comp2[n]', 'met1[n]']))
        self.assertEqual(model.rate_laws.get_one(id='r1_kb-backward').expression.observables,
            [model.observables.get_one(id='obs2')])
        self.assertEqual(model.rate_laws.get_one(id='r1_kb-backward').expression.functions,
            [model.functions.get_one(id='volume_n')])

    def test_unchanged_kb(self):    
        
        kb2 = self.kb.copy()

        self.assertTrue(kb2.is_equal(self.kb))
        
        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel], 
            options={'component': {'InitializeModel': self.set_options([
                        'gen_dna',
                        'gen_transcripts',
                        'gen_protein',
                        'gen_metabolites',
                        'gen_complexes',
                        'gen_distribution_init_concentrations',
                        'gen_observables',
                        'gen_kb_reactions',
                        'gen_kb_rate_laws'])}}).run()
        
        self.assertTrue(kb2.is_equal(self.kb))
