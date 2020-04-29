""" Tests for macromolecular complexation submodel generator

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-08-02
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_model_gen.eukaryote import complexation
from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import wc_model_gen.global_vars as gvar
import os
import scipy.constants
import shutil
import tempfile
import unittest
import wc_kb
import wc_lang
import wc_model_gen.utils as utils


class TestCase(unittest.TestCase):

    def setUp(self):

        # Create KB content
        self.tmp_dirname = tempfile.mkdtemp()
        self.sequence_path = os.path.join(self.tmp_dirname, 'test_seq.fasta')
        with open(self.sequence_path, 'w') as f:
            f.write('>chr1\nGCGTGCGATGATtgatga\n')

        self.kb = wc_kb.KnowledgeBase()
        cell = self.kb.cell = wc_kb.Cell()

        nucleus = cell.compartments.create(id='n')
        mito = cell.compartments.create(id='m')
        lysosome = cell.compartments.create(id='l')
        membrane = cell.compartments.create(id='c_m')
        cytoplasm = cell.compartments.create(id='c')

        chr1 = wc_kb.core.DnaSpeciesType(cell=cell, id='chr1', sequence_path=self.sequence_path)
        gene1 = wc_kb.eukaryote.GeneLocus(cell=cell, id='gene1', polymer=chr1, start=1, end=18)
        
        locus1 = wc_kb.eukaryote.GenericLocus(start=1, end=6)
        transcript1 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, gene=gene1, exons=[locus1])
        prot1 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot1', name='protein1', transcript=transcript1, coding_regions=[locus1])
        prot1_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot1, 
            value='40000.0', value_type=wc_ontology['WC:float'])
        prot1_spec = wc_kb.core.Species(species_type=prot1, compartment=nucleus)
        
        locus2 = wc_kb.eukaryote.GenericLocus(start=4, end=9)
        transcript2 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, gene=gene1, exons=[locus2])
        prot2 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot2', name='protein2', transcript=transcript2, coding_regions=[locus2])
        prot2_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot2, 
            value='20000.0', value_type=wc_ontology['WC:float'])
        prot2_spec = wc_kb.core.Species(species_type=prot2, compartment=nucleus)
        
        locus3 = wc_kb.eukaryote.GenericLocus(start=7, end=12)
        transcript3 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, gene=gene1, exons=[locus3])
        prot3 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot3', name='protein3', transcript=transcript3, coding_regions=[locus3])
        prot3_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot3, 
            value='25000.0', value_type=wc_ontology['WC:float'])
        prot3_spec1 = wc_kb.core.Species(species_type=prot3, compartment=nucleus)
        prot3_spec2 = wc_kb.core.Species(species_type=prot3, compartment=mito)
        prot3_spec3 = wc_kb.core.Species(species_type=prot3, compartment=membrane)

        locus4 = wc_kb.eukaryote.GenericLocus(start=10, end=18)
        transcript4 = wc_kb.eukaryote.TranscriptSpeciesType(cell=cell, gene=gene1, exons=[locus4])
        prot4 = wc_kb.eukaryote.ProteinSpeciesType(cell=cell, id='prot4', name='protein4', transcript=transcript4, coding_regions=[locus4])
        prot4_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=prot4, 
            value='40000.0', value_type=wc_ontology['WC:float'])
        prot4_spec1 = wc_kb.core.Species(species_type=prot4, compartment=nucleus)
        prot4_spec2 = wc_kb.core.Species(species_type=prot4, compartment=cytoplasm)

        trans5 = wc_kb.eukaryote.TranscriptSpeciesType(id='trans5', cell=cell, gene=gene1)
        trans5_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=trans5, 
            value='40000.0', value_type=wc_ontology['WC:float'])
        trans5_spec = wc_kb.core.Species(species_type=trans5, compartment=cytoplasm)

        trans6 = wc_kb.eukaryote.TranscriptSpeciesType(id='trans6', cell=cell, gene=gene1)
        trans6_half_life = wc_kb.core.SpeciesTypeProperty(property='half-life', species_type=trans6, 
            value='36000.0', value_type=wc_ontology['WC:float'])
        trans6_spec = wc_kb.core.Species(species_type=trans6, compartment=cytoplasm)

        met1 = wc_kb.core.MetaboliteSpeciesType(cell=cell, id='met1', name='metabolite1')
        for i in cell.compartments:
            met1_species = wc_kb.core.Species(species_type=met1, compartment=i)
        
        complex1 = wc_kb.core.ComplexSpeciesType(cell=cell, id='complex_1', subunits=[
            wc_kb.core.SpeciesTypeCoefficient(species_type=prot1, coefficient=1),
            wc_kb.core.SpeciesTypeCoefficient(species_type=prot2, coefficient=2),
            wc_kb.core.SpeciesTypeCoefficient(species_type=prot3, coefficient=0),
            ])

        complex2 = wc_kb.core.ComplexSpeciesType(cell=cell, id='complex_2', subunits=[
            wc_kb.core.SpeciesTypeCoefficient(species_type=prot3, coefficient=2),
            wc_kb.core.SpeciesTypeCoefficient(species_type=met1, coefficient=2),
            ])

        complex3 = wc_kb.core.ComplexSpeciesType(cell=cell, id='complex_3', subunits=[
            wc_kb.core.SpeciesTypeCoefficient(species_type=prot4, coefficient=2),
            ])

        complex4 = wc_kb.core.ComplexSpeciesType(cell=cell, id='complex_4', subunits=[
            wc_kb.core.SpeciesTypeCoefficient(species_type=prot4, coefficient=1),
            wc_kb.core.SpeciesTypeCoefficient(species_type=trans5, coefficient=1),
            wc_kb.core.SpeciesTypeCoefficient(species_type=trans6, coefficient=2),
            ])                    

        # Create initial model content
        self.model = model = wc_lang.Model()
        
        model.parameters.create(id='Avogadro', value = scipy.constants.Avogadro,
                                units = unit_registry.parse_units('molecule mol^-1'))

        compartments = {'n': ('nucleus', 5E-14), 'm': ('mitochondria', 2.5E-14), 
            'l': ('lysosome', 2.5E-14), 'c_m': ('membrane', 5E-15), 'c': ('cytoplasm', 1E-13)}
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

        for i in cell.species_types.get(__type=wc_kb.eukaryote.ProteinSpeciesType):
            model_species_type = model.species_types.get_or_create(id=i.id, name=i.name, type=wc_ontology['WC:protein'])
            model_compartment = model.compartments.get_one(id='n')
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.get_or_create(species=model_species, 
                mean=10, units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()

        model_species_type = model.species_types.get_or_create(id='prot3', name='protein3', type=wc_ontology['WC:protein'])
        model_mito = model.compartments.get_one(id='m')
        model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_mito)
        model_species.id = model_species.gen_id()
        conc_model = model.distribution_init_concentrations.create(species=model_species, 
            mean=20, units=unit_registry.parse_units('molecule'))
        conc_model.id = conc_model.gen_id()        
        
        model_membrane = model.compartments.get_one(id='c_m')
        model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_membrane)
        model_species.id = model_species.gen_id()
        conc_model = model.distribution_init_concentrations.create(species=model_species, 
            mean=20, units=unit_registry.parse_units('molecule'))
        conc_model.id = conc_model.gen_id()

        model_species_type = model.species_types.get_or_create(id='prot4', name='protein4', type=wc_ontology['WC:protein'])
        model_cyto = model.compartments.get_one(id='c')
        model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_cyto)
        model_species.id = model_species.gen_id()
        conc_model = model.distribution_init_concentrations.create(species=model_species, 
            mean=20, units=unit_registry.parse_units('molecule'))
        conc_model.id = conc_model.gen_id()

        for i in cell.species_types.get(__type=wc_kb.eukaryote.TranscriptSpeciesType):
            model_species_type = model.species_types.get_or_create(id=i.id, name=i.name, type=wc_ontology['WC:RNA'])
            model_compartment = model.compartments.get_one(id='c')
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()
            conc_model = model.distribution_init_concentrations.get_or_create(species=model_species, 
                mean=25, units=unit_registry.parse_units('molecule'))
            conc_model.id = conc_model.gen_id()

        model_species_type = model.species_types.get_or_create(id='met1', name='metabolite1', type=wc_ontology['WC:metabolite'])
        for compartment in model.compartments:
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=compartment)
            model_species.id = model_species.gen_id()

        for compl in [complex1, complex2, complex3, complex4]:
            model_species_type = model.species_types.create(id=compl.id, type=wc_ontology['WC:pseudo_species'])
            subunit_compartments = [[s.compartment.id for s in sub.species_type.species]
                for sub in compl.subunits]
            if len(subunit_compartments) == 1:
                shared_compartments = set(subunit_compartments[0])
            else:    
                shared_compartments = set([])
                for i in range(len(subunit_compartments)):
                    shared_compartments = (set(subunit_compartments[i])
                        if i==0 else shared_compartments).intersection(
                        set(subunit_compartments[i+1]) if i<(len(subunit_compartments)-1) else shared_compartments)            
            compartment_ids = set(list(shared_compartments))
            for compartment_id in compartment_ids:
                model_compartment = model.compartments.get_one(id=compartment_id)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()

        metabolic_participants = ['amp', 'cmp', 'gmp', 'ump', 'h2o', 'h', 'Ala', 'Cys', 'Asp', 'Glu', 'Selcys']
        metabolic_compartments = ['l', 'm', 'c']
        for i in metabolic_participants:
            for c in metabolic_compartments:
                model_species_type = model.species_types.get_or_create(id=i, type=wc_ontology['WC:metabolite'])            
                model_compartment = model.compartments.get_one(id=c)
                model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
                model_species.id = model_species.gen_id()
            
    def tearDown(self):
        shutil.rmtree(self.tmp_dirname)
        gvar.protein_aa_usage = {}
        gvar.transcript_ntp_usage = {}                     

    def test_methods(self):

        model = self.model

        amino_acid_id_conversion = {
            'A': 'Ala',
            'C': 'Cys',
            'D': 'Asp',
            'E': 'Glu',
            'U': 'Selcys',
            }
        gen = complexation.ComplexationSubmodelGenerator(self.kb, self.model, options={
            'rna_subunit_seq': {'trans5': 'ACCG', 'trans6': 'AU'},
            'amino_acid_id_conversion': amino_acid_id_conversion,
            'cds': False,
            'estimate_initial_state': False,
            'selenoproteome': ['gene1'],                
            })
        gen.run()

        self.assertEqual(gvar.protein_aa_usage['prot4'], 
            {'A': 0, 'C': 0, 'D': 1, 'E': 0, 'U': 1, 'len': 2, '*': 2, 'start_aa': 'D', 'start_codon': 'GAU'})
        self.assertEqual(gvar.transcript_ntp_usage['trans5'], {'A': 1, 'C': 2, 'G': 1, 'U': 0, 'len': 4})       

        self.assertEqual(len(model.reactions), 18)
        self.assertEqual([i.id for i in model.submodels], ['complexation'])
        self.assertEqual(model.submodels.get_one(id='complexation').framework, wc_ontology['WC:next_reaction_method'])

        # Test gen_reaction
        complex1_assembly = model.reactions.get_one(id='complex_1_association_in_n')
        self.assertEqual(complex1_assembly.name, 'Complexation of complex_1 in nucleus')
        self.assertEqual(complex1_assembly.reversible, False)
        self.assertEqual(complex1_assembly.comments, '')
        self.assertEqual([(i.species.id, i.coefficient) for i in complex1_assembly.participants],
            [('prot1[n]', -1), ('prot2[n]', -2), ('prot3[n]', -1), ('complex_1[n]', 1)])
        
        dissociate_prot1 = model.reactions.get_one(id='complex_1_dissociation_in_n_degrade_prot1')
        self.assertEqual(dissociate_prot1.name, 'Dissociation of complex_1 in nucleus and degradation of prot1')
        self.assertEqual(dissociate_prot1.reversible, False)
        self.assertEqual(dissociate_prot1.comments, '')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in dissociate_prot1.participants]),
            sorted([('complex_1[n]', -1), ('h2o[l]', -1), ('Ala[l]', 1), ('Cys[l]', 1), ('prot2[n]', 2), ('prot3[n]', 1)]))
        
        dissociate_prot2 = model.reactions.get_one(id='complex_1_dissociation_in_n_degrade_prot2')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in dissociate_prot2.participants]),
            sorted([('complex_1[n]', -1), ('h2o[l]', -1), ('Cys[l]', 1), ('Asp[l]', 1), ('prot1[n]', 1), ('prot2[n]', 1), ('prot3[n]', 1)]))

        dissociate_prot3 = model.reactions.get_one(id='complex_1_dissociation_in_n_degrade_prot3')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in dissociate_prot3.participants]),
            sorted([('complex_1[n]', -1), ('h2o[l]', -1), ('Asp[l]', 2), ('prot1[n]', 1), ('prot2[n]', 2)]))

        complex2_assembly = model.reactions.get_one(id='complex_2_association_in_m')
        self.assertEqual(complex2_assembly.name, 'Complexation of complex_2 in mitochondria')
        self.assertEqual(complex2_assembly.reversible, False)
        self.assertEqual(complex2_assembly.comments, '')
        self.assertEqual([(i.species.id, i.coefficient) for i in complex2_assembly.participants],
            [('prot3[m]', -2), ('met1[m]', -2), ('complex_2[m]', 1)])

        c2_dissociate_prot3_m = model.reactions.get_one(id='complex_2_dissociation_in_m_degrade_prot3')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in c2_dissociate_prot3_m.participants]),
            sorted([('complex_2[m]', -1), ('h2o[m]', -1), ('Asp[m]', 2),  ('prot3[m]', 1), ('met1[m]', 2)]))

        c2_dissociate_prot3_c_m = model.reactions.get_one(id='complex_2_dissociation_in_c_m_degrade_prot3')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in c2_dissociate_prot3_c_m.participants]),
            sorted([('complex_2[c_m]', -1), ('h2o[l]', -1), ('Asp[l]', 2),  ('prot3[c_m]', 1), ('met1[c_m]', 2)]))

        dissociate_prot4 = model.reactions.get_one(id='complex_3_dissociation_in_n_degrade_prot4')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in dissociate_prot4.participants]),
            sorted([('complex_3[n]', -1), ('h2o[l]', -1), ('Asp[l]', 1), ('Selcys[l]', 1), ('prot4[n]', 1)]))

        complex4_assembly = model.reactions.get_one(id='complex_4_association_in_c')
        self.assertEqual(complex4_assembly.name, 'Complexation of complex_4 in cytoplasm')
        self.assertEqual([(i.species.id, i.coefficient) for i in complex4_assembly.participants],
            [('prot4[c]', -1), ('trans5[c]', -1), ('trans6[c]', -2), ('complex_4[c]', 1)])

        c4_dissociate_prot4_c = model.reactions.get_one(id='complex_4_dissociation_in_c_degrade_prot4')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in c4_dissociate_prot4_c.participants]),
            sorted([('complex_4[c]', -1), ('h2o[l]', -1), ('Asp[l]', 1), ('Selcys[l]', 1), ('trans5[c]', 1), ('trans6[c]', 2)]))

        c4_dissociate_trans5_c = model.reactions.get_one(id='complex_4_dissociation_in_c_degrade_trans5')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in c4_dissociate_trans5_c.participants]),
            sorted([('complex_4[c]', -1), ('h2o[c]', -3), ('amp[c]', 1), ('cmp[c]', 2), ('gmp[c]', 1), ('ump[c]', 0), ('h[c]', 3), 
                ('prot4[c]', 1), ('trans6[c]', 2)]))

        c4_dissociate_trans6_c = model.reactions.get_one(id='complex_4_dissociation_in_c_degrade_trans6')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in c4_dissociate_trans6_c.participants]),
            sorted([('complex_4[c]', -1), ('h2o[c]', -1), ('amp[c]', 1), ('cmp[c]', 0), ('gmp[c]', 0), ('ump[c]', 1), ('h[c]', 1), 
                ('prot4[c]', 1), ('trans5[c]', 1), ('trans6[c]', 1)]))              

        # Test gen_rate_laws
        self.assertEqual(complex1_assembly.rate_laws[0].expression.expression, 
            'k_cat_complex_1_association_in_n * '
            '(prot1[n] / (prot1[n] + K_m_complex_1_association_in_n_prot1 * Avogadro * volume_n)) * '
            '(prot2[n] / (prot2[n] + K_m_complex_1_association_in_n_prot2 * Avogadro * volume_n)) * '
            '(prot3[n] / (prot3[n] + K_m_complex_1_association_in_n_prot3 * Avogadro * volume_n)) * 2 ** 3')
        self.assertEqual(complex1_assembly.rate_laws[0].direction, wc_lang.RateLawDirection.forward)
        self.assertEqual(dissociate_prot1.rate_laws[0].expression.expression, 
            'k_cat_complex_1_dissociation_in_n_degrade_prot1 * complex_1[n]')
        self.assertEqual(dissociate_prot3.rate_laws[0].expression.expression, 
            'k_cat_complex_1_dissociation_in_n_degrade_prot3 * complex_1[n]')
        self.assertEqual(complex4_assembly.rate_laws[0].expression.expression, 
            'k_cat_complex_4_association_in_c * '
            '(prot4[c] / (prot4[c] + K_m_complex_4_association_in_c_prot4 * Avogadro * volume_c)) * '
            '(trans5[c] / (trans5[c] + K_m_complex_4_association_in_c_trans5 * Avogadro * volume_c)) * '
            '(trans6[c] / (trans6[c] + K_m_complex_4_association_in_c_trans6 * Avogadro * volume_c)) * 2 ** 3')
        self.assertEqual(c4_dissociate_trans6_c.rate_laws[0].expression.expression, 
            'k_cat_complex_4_dissociation_in_c_degrade_trans6 * complex_4[c]')

        for law in model.rate_laws:
            self.assertEqual(law.validate(), None)

        # Test calibrate_submodels
        self.assertEqual(model.parameters.get_one(id='K_m_complex_1_association_in_n_prot1').value, 10/scipy.constants.Avogadro/5E-14)
        self.assertEqual(model.parameters.get_one(id='K_m_complex_1_association_in_n_prot3').value, 10/scipy.constants.Avogadro/5E-14)
        self.assertEqual(model.parameters.get_one(id='K_m_complex_1_association_in_n_prot3').comments, 
            'The value was assumed to be 1.0 times the concentration of prot3 in nucleus')        
        self.assertEqual(model.parameters.get_one(id='k_cat_complex_1_association_in_n').value, (1/40000 + 2/20000 + 1/25000)/0.1)
        self.assertEqual(model.parameters.get_one(id='k_cat_complex_1_association_in_n').comments,
            'The value was assigned so that the ratio of effective dissociation constant to '
            'association constant is the same as the specified ratio of free subunit to subunit in complexes '
            'at equilibrium') 
        self.assertEqual(model.parameters.get_one(id='k_cat_complex_1_dissociation_in_n_degrade_prot1').value, 1/40000.)
        self.assertEqual(model.parameters.get_one(id='k_cat_complex_1_dissociation_in_n_degrade_prot2').value, 2/20000.)
        self.assertEqual(model.parameters.get_one(id='k_cat_complex_1_dissociation_in_n_degrade_prot3').value, 1/25000.)
        self.assertEqual(model.parameters.get_one(id='K_m_complex_4_association_in_c_trans5').value, 25/scipy.constants.Avogadro/1E-13)
        self.assertEqual(model.parameters.get_one(id='k_cat_complex_4_dissociation_in_c_degrade_trans6').value, 2/36000.)

        # Test determine_initial_concentration
        gen.options['estimate_initial_state'] = True
        gen.options['greedy_step_size'] = 0.5
        gen.calibrate_submodel()
        self.assertEqual(gen._maximum_possible_amount, {'complex_1[n]': 5.0, 'complex_2[n]': 5.0, 'complex_2[m]': 10.0, 
            'complex_2[c_m]': 10.0, 'complex_3[n]': 5.0, 'complex_3[c]': 10.0, 'complex_4[c]': 12.5})
        print(gen._effective_dissociation_constant)
        self.assertEqual(gen._effective_dissociation_constant, {
            'complex_1[n]': 1/40000 + 2/20000 + 1/25000, 
            'complex_2[n]': 2/25000, 
            'complex_2[m]': 2/25000, 
            'complex_2[c_m]': 2/25000,
            'complex_3[n]': 2/40000,
            'complex_3[c]': 2/40000,
            'complex_4[c]': 1/40000 + 1/40000 + 2/36000,
            })
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-prot1[n]').mean, 6)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-prot2[n]').mean, 2)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-prot3[n]').mean, 1)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-prot3[m]').mean, 2)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-prot3[c_m]').mean, 2)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-complex_1[n]').mean, 4)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-complex_2[n]').mean, 2.5)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-complex_2[m]').mean, 9)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-complex_2[c_m]').mean, 9)
        self.assertEqual(model.distribution_init_concentrations.get_one(id='dist-init-conc-complex_3[n]').mean, 4.5)
        
        self.assertEqual(model.parameters.get_one(id='k_cat_complex_1_association_in_n').value, (1/40000 + 2/20000 + 1/25000) * 4)
        self.assertEqual(model.parameters.get_one(id='k_cat_complex_2_association_in_n').value, 2/25000 * 2.5)
        self.assertEqual(model.parameters.get_one(id='K_m_complex_1_association_in_n_prot1').value, 6/scipy.constants.Avogadro/5E-14)
        self.assertEqual(model.parameters.get_one(id='K_m_complex_1_association_in_n_prot1').comments, 
            'The value was assumed to be 1.0 times the concentration of prot1 in nucleus')
        self.assertEqual(model.parameters.get_one(id='K_m_complex_2_association_in_m_met1').value, 1e-05)
        self.assertEqual(model.parameters.get_one(id='K_m_complex_2_association_in_m_met1').comments, 
            'The value was assigned to 1e-05 because the concentration of met1 in mitochondria was not known')

    def test_global_vars(self):
        gvar.protein_aa_usage = {'prot1': {'A': 4, 'C': 2, 'D': 1, 'E': 0, 'U': 0, 'len': 7, '*': 1, 'start_aa': 'A', 'start_codon': 'GCG'}}
        gvar.transcript_ntp_usage = {
            'trans5': {'A': 1, 'C': 1, 'G': 1, 'U': 1, 'len': 4}, 
            'trans6': {'A': 0, 'C': 1, 'G': 1, 'U': 2, 'len': 4},
            }
        amino_acid_id_conversion = {
            'A': 'Ala',
            'C': 'Cys',
            'D': 'Asp',
            'E': 'Glu',
            'U': 'Selcys',
            }
        gen = complexation.ComplexationSubmodelGenerator(self.kb, self.model, options={
            'amino_acid_id_conversion': amino_acid_id_conversion,
            'cds': False,
            'estimate_steady_state': False,
            'selenoproteome': ['gene1'],                
            })
        gen.run()   

        dissociate_prot1 = self.model.reactions.get_one(id='complex_1_dissociation_in_n_degrade_prot1')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in dissociate_prot1.participants]),
            sorted([('complex_1[n]', -1), ('h2o[l]', -6), ('Ala[l]', 4), ('Cys[l]', 2), ('Asp[l]', 1), ('prot2[n]', 2), ('prot3[n]', 1)]))
        
        c4_dissociate_trans6_c = self.model.reactions.get_one(id='complex_4_dissociation_in_c_degrade_trans6')
        self.assertEqual(sorted([(i.species.id, i.coefficient) for i in c4_dissociate_trans6_c.participants]),
            sorted([('complex_4[c]', -1), ('h2o[c]', -3), ('amp[c]', 0), ('cmp[c]', 1), ('gmp[c]', 1), ('ump[c]', 2), ('h[c]', 3), 
                ('prot4[c]', 1), ('trans5[c]', 1), ('trans6[c]', 1)]))
