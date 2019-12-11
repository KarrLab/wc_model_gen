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
import wc_model_gen.global_vars as gvar
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
                        'gen_environment': False,
                        }

        for i in options:
            del option_dict[i]

        return option_dict

    def setUp(self):

        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>chr1\nTTTatgaARGTNCTCATHAAYAARAAYGARCTCTAGTTTAT\n'
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
        cytoplasm = cell.compartments.create(id='c', name='cytoplasm', volumetric_fraction=0.5)
        mito = cell.compartments.create(id='m', name='mitochondrion', volumetric_fraction=0.2)
        extra = cell.compartments.create(id='e', name='extracellular space')

        chr1 = wc_kb.core.DnaSpeciesType(cell=cell, id='chr1', name='chromosome 1', ploidy=2,
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene1 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene1', polymer=chr1, start=1, end=36)
        exon1 = wc_kb.eukaryote.GenericLocus(start=4, end=36)
        transcript1 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans1',
            name='transcript1', gene=gene1, exons=[exon1])
        transcript1_spec = wc_kb.core.Species(species_type=transcript1, compartment=cytoplasm)
        transcript1_conc = wc_kb.core.Concentration(cell=cell, species=transcript1_spec, value=0.02)
        transcript1_conc.id = transcript1_conc.serialize()
        cds1 = wc_kb.eukaryote.GenericLocus(start=4, end=36)        
        prot1 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot1', name='protein1', 
            transcript=transcript1, coding_regions=[cds1])
        prot1_spec = wc_kb.core.Species(species_type=prot1, compartment=nucleus)
        prot1_conc = wc_kb.core.Concentration(cell=cell, species=prot1_spec, value=0.03)
        prot1_conc.id = prot1_conc.serialize()

        chrX = wc_kb.core.DnaSpeciesType(cell=cell, id='chrX', name='chromosome X', ploidy=1, 
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene2 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene2', polymer=chrX, start=1, end=4)
        exon2 = wc_kb.eukaryote.GenericLocus(start=1, end=4)
        transcript2 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans2',
            name='transcript2', gene=gene2, exons=[exon2])
        transcript2_spec = wc_kb.core.Species(species_type=transcript2, compartment=cytoplasm)
        transcript2_conc = wc_kb.core.Concentration(cell=cell, species=transcript2_spec, value=0.01)
        transcript2_conc.id = transcript2_conc.serialize()         

        chrM = wc_kb.core.DnaSpeciesType(cell=cell, id='chrM', name='mitochondrial chromosome', ploidy=150,
            sequence_path=self.sequence_path, circular=False, double_stranded=False)
        gene3 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene3', polymer=chrM, start=1, end=33)
        exon3 = wc_kb.eukaryote.GenericLocus(start=1, end=30)
        transcript3 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, id='trans3',
            name='transcript3', gene=gene3, exons=[exon3])
        transcript3_spec = wc_kb.core.Species(species_type=transcript3, compartment=mito)
        transcript3_conc = wc_kb.core.Concentration(cell=cell, species=transcript3_spec, value=0.05)
        transcript3_conc.id = transcript3_conc.serialize()
        cds3 = wc_kb.eukaryote.GenericLocus(start=1, end=30)        
        prot3 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot3', name='protein3', 
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

        electron = wc_kb.core.MetaboliteSpeciesType(cell=cell, id='el', name='electron')
        proton = wc_kb.core.MetaboliteSpeciesType(cell=cell, id='h', name='proton')
        pmf = wc_kb.core.MetaboliteSpeciesType(cell=cell, id='PMF', name='proton motive force')
        h2 = wc_kb.core.MetaboliteSpeciesType(cell=cell, id='h2', name='dihydrogen')

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

        species_type_coeff5 = wc_kb.core.SpeciesTypeCoefficient(species_type=prot1, coefficient=0)
        complex3 = wc_kb.core.ComplexSpeciesType(cell=cell, id='comp3', name='complex3',
            subunits=[species_type_coeff5])

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

        met2 = wc_kb.core.MetaboliteSpeciesType(cell=cell, id='met2', name='metabolite2')
        met2_formula = wc_kb.core.SpeciesTypeProperty(property='empirical_formula', species_type=met2, 
            value='H3', value_type=wc_ontology['WC:string'])
        met2_formula.id = met2_formula.gen_id()
        met2_charge = wc_kb.core.SpeciesTypeProperty(property='charge', species_type=met2, 
            value=3, value_type=wc_ontology['WC:integer'])
        met2_charge.id = met2_charge.gen_id()
        proton_species = wc_kb.core.Species(species_type=proton, compartment=nucleus)
        met2_species = wc_kb.core.Species(species_type=met2, compartment=nucleus)
        proton_substrate = wc_kb.core.SpeciesCoefficient(species=proton_species, coefficient=-1)
        proton_product = wc_kb.core.SpeciesCoefficient(species=proton_species, coefficient=1)
        met2_substrate = wc_kb.core.SpeciesCoefficient(species=met2_species, coefficient=-1)
        met2_product = wc_kb.core.SpeciesCoefficient(species=met2_species, coefficient=1)        
        reaction2 = wc_kb.core.Reaction(cell=cell, id='r2', name='reaction2',
            participants=[proton_substrate, met2_product], reversible=True)
        reaction3 = wc_kb.core.Reaction(cell=cell, id='r3', name='reaction3',
            participants=[met2_substrate, proton_product], reversible=True)
        reaction4 = wc_kb.core.Reaction(cell=cell, id='r4', name='reaction4',
            participants=[met2_substrate], reversible=True) 
        reaction5 = wc_kb.core.Reaction(cell=cell, id='r5', name='reaction5',
            participants=[participant1], reversible=True)
        reaction6 = wc_kb.core.Reaction(cell=cell, id='r6', name='reaction6',
            participants=[participant2], reversible=True)                       

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
        gvar.protein_aa_usage = {} 

    def test_gen_taxon_compartments_parameters(self):

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([])}}).run()

        self.assertEqual(model.taxon.id, 'taxon')
        self.assertEqual(model.taxon.name, 'Homo sapiens')
        self.assertEqual(model.taxon.rank, wc_lang.TaxonRank.species)

        self.assertEqual(model.parameters.get_one(id='Avogadro').value, scipy.constants.Avogadro)
        self.assertEqual(model.parameters.get_one(id='Avogadro').type, None)
        self.assertEqual(model.parameters.get_one(id='Avogadro').units, unit_registry.parse_units('molecule mol^-1'))

        self.assertEqual(model.compartments.get_one(id='n').name, 'nucleus')
        self.assertEqual(model.compartments.get_one(id='n').biological_type, wc_ontology['WC:cellular_compartment'])
        self.assertAlmostEqual(model.compartments.get_one(id='n').init_volume.mean, 5199.999998548484)
        self.assertAlmostEqual(model.compartments.get_one(id='n_m').init_volume.mean, 1.4515160909356243e-06)
        self.assertEqual(model.compartments.get_one(id='n_m').physical_type, wc_ontology['WC:membrane_compartment'])
        self.assertEqual(model.compartments.get_one(id='e').init_volume.mean, 1.0)
        self.assertEqual(model.compartments.get_one(id='e').biological_type, wc_ontology['WC:extracellular_compartment'])

        self.assertEqual(model.parameters.get_one(id='density_n').value, 1040.)
        self.assertEqual(model.parameters.get_one(id='density_n_m').value, 1160.)
        self.assertEqual(model.parameters.get_one(id='density_e').value, 1000.)
        self.assertEqual(model.parameters.get_one(id='density_n').units, unit_registry.parse_units('g l^-1'))

        self.assertEqual(model.functions.get_one(id='volume_n').units, unit_registry.parse_units('l'))
        self.assertEqual(wc_lang.Function.expression.serialize(
            model.functions.get_one(id='volume_n').expression), 'n / density_n')
        self.assertEqual(wc_lang.Function.expression.serialize(
            model.functions.get_one(id='volume_n_m').expression), 'n_m / density_n_m')
        self.assertEqual(wc_lang.Function.expression.serialize(
            model.functions.get_one(id='volume_e').expression), 'e / density_e')
        
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').value, 0.2)
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').units, unit_registry.parse_units('molecule^-1 s^-1'))

        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').type, wc_ontology['WC:k_cat'])
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').identifiers[0].namespace, 'Sabio-RK')
        self.assertEqual(model.parameters.get_one(id='k_cat_r1_forward').identifiers[0].id, '1234')
        self.assertEqual(model.parameters.get_one(id='K_m_r1_backward_met1_n').value, 0.3)
        self.assertEqual(model.parameters.get_one(id='K_m_r1_backward_met1_n').units, unit_registry.parse_units('M'))
        self.assertEqual(model.parameters.get_one(id='K_m_r1_backward_met1_n').type, wc_ontology['WC:K_m'])

        self.assertEqual(model.parameters.get_one(id='cell_volume').references[0].id, 'ref_1')
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
        gvar.transcript_ntp_usage = {}
        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_transcripts'])}}).run()

        transcript1_model = model.species_types.get_one(id='trans1')
        transcript2_model = model.species_types.get_one(id='trans2')
        transcript3_model = model.species_types.get_one(id='trans3')

        self.assertEqual(gvar.transcript_ntp_usage['trans2'], {'A': 1, 'U': 1, 'G': 1, 'C': 1, 'len': 4})

        self.assertEqual(transcript1_model.name, 'transcript1')
        self.assertEqual(transcript3_model.type, wc_ontology['WC:RNA'])
        self.assertEqual(all(i.compartment.id=='c' for i in model.species.get(species_type=transcript1_model)), True)
        self.assertEqual(all(i.compartment.id=='c' for i in model.species.get(species_type=transcript2_model)), True)
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

        self.assertEqual(gvar.protein_aa_usage['prot1'], {
                'len': 10,
                '*': 0,  # Symbol used in Bio.Seq.Seq when cds is set to False  
                'A': 0,  # Ala: Alanine (C3 H7 N O2)
                'R': 0,  # Arg: Arginine (C6 H14 N4 O2)
                'N': 2,  # Asn: Asparagine (C4 H8 N2 O3)
                'D': 0,  # Asp: Aspartic acid (C4 H7 N O4)
                'C': 0,  # Cys: Cysteine (C3 H7 N O2 S)
                'Q': 0,  # Gln: Glutamine (C5 H10 N2 O3)
                'E': 1,  # Glu: Glutamic acid (C5 H9 N O4)
                'G': 0,  # Gly: Glycine (C2 H5 N O2)
                'H': 0,  # His: Histidine (C6 H9 N3 O2)
                'I': 1,  # Ile: Isoleucine (C6 H13 N O2)
                'L': 2,  # Leu: Leucine (C6 H13 N O2)
                'K': 2,  # Lys: Lysine (C6 H14 N2 O2)
                'M': 1,  # Met: Methionine (C5 H11 N O2 S)
                'F': 0,  # Phe: Phenylalanine (C9 H11 N O2)
                'P': 0,  # Pro: Proline (C5 H9 N O2)
                'S': 0,  # Ser: Serine (C3 H7 N O3)
                'T': 0,  # Thr: Threonine (C4 H9 N O3)
                'W': 0,  # Trp: Tryptophan (C11 H12 N2 O2)
                'Y': 0,  # Tyr: Tyrosine (C9 H11 N O3)
                'V': 1,  # Val: Valine (C5 H11 N O2)
                'U': 0,  # Selcys: Selenocysteine (C3 H7 N O2 Se)
            })

        self.assertEqual(prot1_model.name, 'protein1')
        self.assertEqual(prot3_model.type, wc_ontology['WC:protein'])
        self.assertEqual(all(i.compartment.id=='n' for i in model.species.get(species_type=prot1_model)), True)
        self.assertEqual(all(i.compartment.id=='m' for i in model.species.get(species_type=prot3_model)), True)
        self.assertEqual(prot1_model.structure.empirical_formula, chem.EmpiricalFormula('C53H96N14O15S1'))
        self.assertAlmostEqual(prot1_model.structure.molecular_weight, 1201.49, delta=0.3)
        self.assertEqual(prot1_model.structure.charge, 1)
        self.assertEqual(prot1_model.comments, '')
        self.assertEqual(prot3_model.comments, 'Testing')

        # Test the use of determine_protein_structure_from_aa
        amino_acid_id_conversion = {'N': 'met_N', 'E': 'met_E', 'I': 'met_I', 'L': 'met_L', 'K':'met_K', 'M':'met_M', 'V':'met_V'}
        for i in amino_acid_id_conversion.values():
            kb_met = self.kb.cell.species_types.create(__type=wc_kb.core.MetaboliteSpeciesType, id=i)
            kb_met.properties.create(property='empirical_formula', value='C6H12O6N', value_type=wc_ontology['WC:string'])
            kb_met.properties.create(property='charge', value='-1', value_type=wc_ontology['WC:integer'])
        
        option = self.set_options(['gen_metabolites', 'gen_protein'])
        option.update({'amino_acid_id_conversion': amino_acid_id_conversion})
        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': 
                option}}).run()
        
        prot1_model = model.species_types.get_one(id='prot1')
        self.assertEqual(prot1_model.structure.empirical_formula, chem.EmpiricalFormula('C60H102N10O51'))
        self.assertAlmostEqual(prot1_model.structure.molecular_weight, 1779.4950000000001, delta=0.3)
        self.assertEqual(prot1_model.structure.charge, -10)           

    def test_gen_metabolites(self):

        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options(['gen_metabolites'])}}).run()

        met1_model = model.species_types.get_one(id='met1')

        self.assertEqual(met1_model.name, 'metabolite1')
        self.assertEqual(met1_model.type, wc_ontology['WC:metabolite'])
        self.assertEqual(met1_model.structure.value, 'NC1=C2N=CN([C@@H]3O[C@H](COP([O-])([O-])=O)[C@@H](O)[C@H]3O)C2=NC=N1')
        self.assertEqual(met1_model.structure.format, wc_lang.ChemicalStructureFormat.SMILES)
        self.assertEqual(set([i.compartment.id for i in model.species.get(species_type=met1_model)]), set(['n', 'e']))
        self.assertEqual(met1_model.structure.empirical_formula, chem.EmpiricalFormula('C10H12N5O7P'))
        self.assertAlmostEqual(met1_model.structure.molecular_weight, 345.20776199799997, places=4)
        self.assertEqual(met1_model.structure.charge, -2)
        self.assertEqual(met1_model.comments, '')

        electron_model = model.species_types.get_one(id='el')

        self.assertEqual(electron_model.name, 'electron')
        self.assertEqual(electron_model.type, wc_ontology['WC:metabolite'])
        self.assertEqual(electron_model.structure.value, '[*-]')
        self.assertEqual(electron_model.structure.format, wc_lang.ChemicalStructureFormat.SMILES)
        self.assertEqual(electron_model.structure.empirical_formula, chem.EmpiricalFormula())
        self.assertEqual(electron_model.structure.molecular_weight, 0.)
        self.assertEqual(electron_model.structure.charge, -1)

        proton_model = model.species_types.get_one(id='h')

        self.assertEqual(proton_model.name, 'proton')
        self.assertEqual(proton_model.type, wc_ontology['WC:metabolite'])
        self.assertEqual(proton_model.structure.value, '[1H+]')
        self.assertEqual(proton_model.structure.format, wc_lang.ChemicalStructureFormat.SMILES)
        self.assertEqual(proton_model.structure.empirical_formula, chem.EmpiricalFormula('H'))
        self.assertEqual(proton_model.structure.molecular_weight, 1.008)
        self.assertEqual(proton_model.structure.charge, 1)

        pmf_model = model.species_types.get_one(id='PMF')

        self.assertEqual(pmf_model.name, 'proton motive force')
        self.assertEqual(pmf_model.type, wc_ontology['WC:metabolite'])
        self.assertEqual(pmf_model.structure.value, '[1H+]')
        self.assertEqual(pmf_model.structure.format, wc_lang.ChemicalStructureFormat.SMILES)
        self.assertEqual(pmf_model.structure.empirical_formula, chem.EmpiricalFormula('H'))
        self.assertEqual(pmf_model.structure.molecular_weight, 1.008)
        self.assertEqual(pmf_model.structure.charge, 1)

        h2_model = model.species_types.get_one(id='h2')

        self.assertEqual(h2_model.name, 'dihydrogen')
        self.assertEqual(h2_model.type, wc_ontology['WC:metabolite'])
        self.assertEqual(h2_model.structure.value, '[H][H]')
        self.assertEqual(h2_model.structure.format, wc_lang.ChemicalStructureFormat.SMILES)
        self.assertEqual(h2_model.structure.empirical_formula, chem.EmpiricalFormula('H2'))
        self.assertEqual(h2_model.structure.molecular_weight, 2.016)
        self.assertEqual(h2_model.structure.charge, 0)
        
    def test_gen_complexes(self):

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

        comp3_model = model.species_types.get_one(id='comp3')
        self.assertEqual(set([i.compartment.id for i in model.species.get(species_type=comp3_model)]), set(['n']))
        self.assertEqual(comp3_model.structure.empirical_formula, chem.EmpiricalFormula('C53H96N14O15S1'))
        self.assertAlmostEqual(comp3_model.structure.molecular_weight, 1201.49, delta=0.3)
        self.assertEqual(comp3_model.structure.charge, 1)

        # Test the use of determine_protein_structure_from_aa
        amino_acid_id_conversion = {'N': 'met_N', 'E': 'met_E', 'I': 'met_I', 'L': 'met_L', 'K':'met_K', 'M':'met_M', 'V':'met_V'}
        for i in amino_acid_id_conversion.values():
            kb_met = self.kb.cell.species_types.create(__type=wc_kb.core.MetaboliteSpeciesType, id=i)
            kb_met.properties.create(property='empirical_formula', value='C6H12O6N', value_type=wc_ontology['WC:string'])
            kb_met.properties.create(property='charge', value='-1', value_type=wc_ontology['WC:integer'])
        
        option = self.set_options(['gen_metabolites', 'gen_complexes'])
        option.update({'amino_acid_id_conversion': amino_acid_id_conversion})
        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': 
                option}}).run()
        
        comp3_model = model.species_types.get_one(id='comp3')
        self.assertEqual(comp3_model.structure.empirical_formula, chem.EmpiricalFormula('C60H102N10O51'))
        self.assertAlmostEqual(comp3_model.structure.molecular_weight, 1779.4950000000001, delta=0.3)
        self.assertEqual(comp3_model.structure.charge, -10)

        option = self.set_options(['gen_metabolites', 'gen_protein', 'gen_complexes'])
        option.update({'amino_acid_id_conversion': amino_acid_id_conversion})
        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': 
                option}}).run()
        
        comp3_model = model.species_types.get_one(id='comp3')
        self.assertEqual(comp3_model.structure.empirical_formula, chem.EmpiricalFormula('C60H102N10O51'))
        self.assertAlmostEqual(comp3_model.structure.molecular_weight, 1779.4950000000001, delta=0.3)
        self.assertEqual(comp3_model.structure.charge, -10)

    def test_gen_distribution_init_concentrations(self):

        test_conc = self.kb.cell.concentrations.get_one(value=0.5)
        test_conc.comments = 'Testing'
        test_conc.references.append(wc_kb.core.Reference(id='ref10', title='Title1', authors='Author1', 
            journal='Journal1', volume='1', issue='1', pages='20', year=1999, comments='xyz'))
        test_conc.identifiers.append(wc_kb.core.Identifier(namespace='ECMDB', id='12345'))

        ref1 = wc_lang.Reference(id='ref10', title='Title1', author='Author1', 
            publication='Journal1', volume='1', issue='1', pages='20', year=1999, comments='xyz')
        ref2 = wc_lang.Reference(title='Title2', author='Author2', 
            publication='Journal2', volume='2', issue='2', pages='200', year=1998)      
        options = self.set_options(['gen_dna', 'gen_distribution_init_concentrations'])
        options['media'] = {'met1': (0.75, [ref1, ref2], 'Comments for extracellular concentration')}
        model = core.EukaryoteModelGenerator(self.kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': options}}).run()

        met1_nucleus = model.distribution_init_concentrations.get_one(id='dist-init-conc-met1[n]')

        self.assertEqual(met1_nucleus.species, model.species.get_one(id='met1[n]'))
        self.assertEqual(met1_nucleus.mean, 0.5*scipy.constants.Avogadro*(0.5*10400.-4.836E-09*(0.5*10400.)**(2/3)))
        self.assertEqual(met1_nucleus.units, unit_registry.parse_units('molecule'))
        self.assertEqual(met1_nucleus.comments, 'Testing')
        self.assertEqual(met1_nucleus.references[0].id, 'ref_2')
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

        met1_extra = model.distribution_init_concentrations.get_one(id='dist-init-conc-met1[e]')

        self.assertEqual(met1_extra.mean, 0.75*scipy.constants.Avogadro*1.)
        self.assertEqual(met1_extra.units, unit_registry.parse_units('molecule'))
        self.assertEqual(met1_extra.comments, 'Comments for extracellular concentration')
        self.assertEqual(sorted([i.title for i in met1_extra.references]), sorted(['Title1', 'Title2']))
        self.assertEqual([i.model for i in met1_extra.references], [model, model])

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

        kb = self.kb
        cell = kb.cell

        model = core.EukaryoteModelGenerator(kb,
            component_generators=[initialize_model.InitializeModel],
            options={'component': {'InitializeModel': self.set_options([
                'gen_protein', 'gen_metabolites', 'gen_complexes', 'gen_kb_reactions'])}}).run()

        self.assertEqual(len(model.reactions), 5)
        self.assertEqual(len(model.submodels), 1)

        self.assertEqual(model.reactions.get_one(id='r1_kb').name, 'reaction1')
        self.assertEqual(model.reactions.get_one(id='r1_kb').submodel.id, 'metabolism')
        self.assertEqual(model.reactions.get_one(id='r1_kb').submodel.framework, wc_ontology['WC:dynamic_flux_balance_analysis'])
        self.assertEqual(model.reactions.get_one(id='r1_kb').submodel.dfba_obj.id, 'dfba-obj-metabolism')
        self.assertEqual(model.reactions.get_one(id='r1_kb').reversible, True)
        self.assertEqual(model.reactions.get_one(id='r1_kb').comments, '')
        self.assertEqual(model.reactions.get_one(id='r6_kb'), None)

        attr = wc_lang.core.ReactionParticipantAttribute()
        self.assertEqual(attr.serialize(model.reactions.get_one(id='r1_kb').participants), 'met1[e] ==> met1[n]')            
        self.assertEqual(attr.serialize(model.reactions.get_one(id='r2_kb').participants), '[n]: (3) h ==> met2')            
        self.assertEqual(attr.serialize(model.reactions.get_one(id='r3_kb').participants), '[n]: met2 ==> (3) h')
        self.assertEqual(attr.serialize(model.reactions.get_one(id='r4_kb').participants), '[n]: met2 ==> (3) h')
        self.assertEqual(attr.serialize(model.reactions.get_one(id='r5_kb').participants), '[n]:  ==> met1')
        
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

    def test_gen_environment(self):

        environment = {'id': 'env', 'name': 'test_environment', 'temperature': 37., 'comments': ''}
        test_option = self.set_options(['gen_environment'])
        test_option['environment'] = environment

        model = core.EukaryoteModelGenerator(self.kb, 
            component_generators=[initialize_model.InitializeModel], 
            options={'component': {'InitializeModel': test_option
                }}).run()

        self.assertEqual(model.env.id, 'env')
        self.assertEqual(model.env.name, 'test_environment')
        self.assertEqual(model.env.temp, 37.)
        self.assertEqual(model.env.temp_units, unit_registry.parse_units('celsius'))
        self.assertEqual(model.env.comments, '')

    def test_structure_to_smiles_and_props(self):

        model = wc_lang.Model()
        test_instance = initialize_model.InitializeModel(self.kb, model)
        ph = 7.4        
        structure = 'InChI=1S/C10H14N5O7P/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(22-10)1-21-23(18,19)20' +\
                    '/h2-4,6-7,10,16-17H,1H2,(H2,11,12,13)(H2,18,19,20)/p-2/t4-,6-,7-,10-/m1/s1'
        smiles, formula, charge, mol_wt = test_instance.structure_to_smiles_and_props(structure, ph)
        self.assertEqual(smiles, 'NC1=C2N=CN([C@@H]3O[C@H](COP([O-])([O-])=O)[C@@H](O)[C@H]3O)C2=NC=N1')
        self.assertEqual(formula, chem.EmpiricalFormula('C10H12N5O7P'))
        self.assertEqual(charge, -2)
        self.assertAlmostEqual(mol_wt, 345.20776199799997, places=4)

        structure = '[H]OC([H])([H])[C@@]1([H])O[C@]([H])(n2c([H])[nH]c3c([H])c(c(c([H])c23)C([H])([H])' +\
                    '[H])C([H])([H])[H])[C@]([H])(O[H])[C@]1([H])OP([O-])(=O)O[C@]([H])(C([H])([H])[H])C([H])' +\
                    '([H])N([H])C(=O)C([H])([H])C([H])([H])[C@]1(\C2=C(\C3=N\C(=C([H])/C4=NC(=C(C5=N[C@@](C([H])' +\
                    '([H])[H])([C@]([H])(N2[Co+][O+]([H])[H])[C@]1([H])C([H])([H])C(=O)N([H])[H])[C@@](C([H])([H])' +\
                    '[H])(C([H])([H])C(=O)N([H])[H])[C@]5([H])C([H])([H])C([H])([H])C(=O)N([H])[H])C([H])([H])[H])' +\
                    '[C@@](C([H])([H])[H])(C([H])([H])C(=O)N([H])[H])[C@]4([H])C([H])([H])C([H])([H])C(=O)N([H])[H])\C' +\
                    '(C([H])([H])[H])(C([H])([H])[H])[C@]3([H])C([H])([H])C([H])([H])C(=O)N([H])[H])C([H])([H])[H])C([H])([H])[H]'
        smiles, formula, charge, mol_wt = test_instance.structure_to_smiles_and_props(structure, ph)
        self.assertEqual(smiles, 'C[C@H](CNC(=O)CC[C@]1(C)[C@@H](CC(N)=O)[C@H]2N([Co+]O)\\C1=C(C)/C1=[NH+]/C(=C\\C3=[NH+]C'
            '(=C(C)C4=N[C@]2(C)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)/C(C)(C)[C@@H]1CCC(N)=O)OP'
            '([O-])(=O)O[C@@H]1[C@@H](CO)O[C@@H]([C@@H]1O)n1c[nH]c2cc(C)c(C)cc12')
        self.assertEqual(formula, chem.EmpiricalFormula('C62H99CoN13O15P'))
        self.assertEqual(charge, 2)
        self.assertAlmostEqual(mol_wt, 1356.456955998, places=4)

    def test_determine_protein_structure_from_aa(self):
        
        model = wc_lang.Model()        
        test_instance = initialize_model.InitializeModel(self.kb, model)
        test_instance.clean_and_validate_options()
        count = {'A': 2, 'R': 3}
        
        formula, weight, charge, determined = test_instance.determine_protein_structure_from_aa('prot1', count)
        self.assertEqual((formula, weight, charge, determined), (chem.EmpiricalFormula(), 0, 0, False))
        
        test_instance.options['amino_acid_id_conversion'] = {'A1': 'a1', 'R1': 'r1'}
        formula, weight, charge, determined = test_instance.determine_protein_structure_from_aa('prot1', count)
        self.assertEqual((formula, weight, charge, determined), (chem.EmpiricalFormula(), 0, 0, False))

        amino_acid = {'a': ('C10H20O10', 14, 3), 'r': ('C20H70O20', 31, 0)}
        for k, v in amino_acid.items():
            model.species_types.create(id=k, structure = wc_lang.ChemicalStructure(
                empirical_formula = chem.EmpiricalFormula(v[0]),
                molecular_weight = v[1],
                charge = v[2],
                )
            )                
        test_instance.options['amino_acid_id_conversion'] = {'A': 'a', 'R': 'r'}
        check_formula = chem.EmpiricalFormula('C10H20O10')*2 + chem.EmpiricalFormula('C20H70O20')*3 - chem.EmpiricalFormula('H8O4')
        check_weight = 2*14 + 3*31 - 8*mendeleev.element('H').atomic_weight - 4*mendeleev.element('O').atomic_weight 
        formula, weight, charge, determined = test_instance.determine_protein_structure_from_aa('prot1', count)
        self.assertEqual((formula, weight, charge, determined), (check_formula, check_weight, 2*3+3*0, True))

        prot1 = model.species_types.create(id='prot1', structure=wc_lang.ChemicalStructure())
        self.assertEqual(prot1.structure.empirical_formula, None)
                
        formula, weight, charge, determined = test_instance.determine_protein_structure_from_aa('prot1', count)
        self.assertEqual((formula, weight, charge, determined), (check_formula, check_weight, 2*3+3*0, True))
        self.assertEqual(prot1.structure.empirical_formula, check_formula)
        self.assertEqual(prot1.structure.molecular_weight, check_weight)
        self.assertEqual(prot1.structure.charge, 2*3+3*0)

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
