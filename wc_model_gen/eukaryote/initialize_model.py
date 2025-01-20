""" Initialize the construction of wc_lang-encoded models from wc_kb-encoded knowledge base.

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-09
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.chem import EmpiricalFormula, OpenBabelUtils
from wc_utils.util.chem.marvin import get_major_micro_species
from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
from wc_utils.util import chem
import wc_model_gen.global_vars as gvar
import math
import mendeleev
import numpy
import openbabel
import scipy.constants
import taxoniq
import wc_kb
import wc_lang
import wc_model_gen


class InitializeModel(wc_model_gen.ModelComponentGenerator):
    """ Initialize model from knowledge base 
    
    Options:

    * culture_volume (:obj:`float`): volume of cell culture; default is 1.0 liter
    * cell_density(:obj:`float`): cell density; default is 1040 g/liter
    * membrane_density (:obj:`float`): membrane density; default is 1160 g/liter
    * cds (:obj:`bool`): True indicates mRNA sequence is a complete CDS; default is True
    * amino_acid_id_conversion (:obj:`dict`): a dictionary with amino acid standard ids
        as keys and amino acid metabolite ids as values
    * selenoproteome (:obj:`list`): list of IDs of genes that translate into 
        selenoproteins, default is an empty list  
    * environment (:obj:`dict`): dictionary with details for generating cell environment in the model 
    * ph (:obj:`float`): pH at which species will be protonated and reactions will be balanced; default is 7.4
    * media (:obj:`dict`): a dictionary with species type ids as keys and tuples of concentration (M) in the 
        media (extracellular space), `list` of `wc_lang.Reference`, and comments as values
    * rna_input_seq (:obj:`dict`, optional): a dictionary with RNA ids as keys and sequence strings as values
    * smiles_input (:obj:`dict`, optional): a dictionary with metabolite ids as keys and smiles strings as values     
    * check_reaction (:obj:`bool`): if True, reactions will be checked and corrected for proton and charge balance;
        default is True
    * gen_dna (:obj:`bool`): if True, DNA species types and species will be generated; 
        default is True
    * gen_transcripts (:obj:`bool`): if True, transcript species types and species will be generated; 
        default is True
    * gen_protein (:obj:`bool`): if True, protein species types and species will be generated; 
        default is True
    * gen_metabolites (:obj:`bool`): if True, metabolite species types and species will be generated; 
        default is True
    * gen_complexes (:obj:`bool`): if True, macromolecular complex species types and species will be generated; 
        default is True
    * gen_distribution_init_concentration (:obj:`bool`): if True, initial concentration of species will be generated; 
        default is True                    
    * gen_observables (:obj:`bool`): if True, observables will be generated; default is True    
    * gen_kb_reactions (:obj:`bool`): if True, reactions will be generated; default is True
    * gen_dfba_objective (:obj:`bool`): if True, a dfba objective function will be created; default is False 
    * gen_kb_rate_laws (:obj:`bool`): if True, rate laws will be generated; default is True
    * gen_environment (:obj:`bool`): if True, cell environment will be generated; default is True    
    """

    def run(self):
        """ Run all the components for initializing model from knowledge base """
        self.clean_and_validate_options()
        options = self.options

        print('Initialization is starting...')

        self.gen_taxon()
        self.gen_compartments()
        self.gen_parameters()
        self.global_vars_from_input()

        print('Taxon, compartments, and parameters have been initialized')

        if options['gen_metabolites']:
            self.gen_metabolites()
            print('All metabolite species types and species have been initialized')

        if options['gen_dna']:
            self.gen_dna()
            print('All DNA species types and species have been initialized')    

        if options['gen_transcripts']:
            self.gen_transcripts()
            print('All transcript species types and species have been initialized')    

        if options['gen_protein']:
            self.gen_protein()
            print('All protein species types and species have been initialized')                

        if options['gen_complexes']:
            self.gen_complexes()
            print('All complex species types and species have been initialized')    

        if options['gen_distribution_init_concentrations']:
            self.gen_distribution_init_concentrations()
            print('Concentrations of species have been initialized')    

        if options['gen_observables']:
            self.gen_observables()
            print('Model observables have been initialized')    

        if options['gen_kb_reactions']:
            self.gen_kb_reactions()
            print('Reactions in knowledge base have been initialized')    

        if options['gen_kb_rate_laws']:
            self.gen_kb_rate_laws()
            print('Rate laws in knowledge base have been initialized')    

        if options['gen_environment']:
            self.gen_environment()
            
        print('Model generator has been initialized')        

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        culture_volume = options.get('culture_volume', 1.)
        assert(isinstance(culture_volume, float))
        options['culture_volume'] = culture_volume

        cell_density = options.get('cell_density', 1040.)
        assert(isinstance(cell_density, float))
        options['cell_density'] = cell_density

        membrane_density = options.get('membrane_density', 1160.)
        assert(isinstance(membrane_density, float))
        options['membrane_density'] = membrane_density

        cds = options.get('cds', True)
        assert(isinstance(cds, bool))
        options['cds'] = cds

        amino_acid_id_conversion = options.get('amino_acid_id_conversion', {})
        assert(isinstance(amino_acid_id_conversion, dict))
        options['amino_acid_id_conversion'] = amino_acid_id_conversion

        selenoproteome = options.get('selenoproteome', [])
        assert(isinstance(selenoproteome, list))
        options['selenoproteome'] = selenoproteome

        environment = options.get('environment', {})
        assert(isinstance(environment, dict))
        options['environment'] = environment

        ph = options.get('ph', 7.4)
        assert(isinstance(ph, float))
        options['ph'] = ph

        media = options.get('media', {})
        assert(isinstance(media, dict))
        options['media'] = media

        rna_input_seq = options.get('rna_input_seq', {})
        assert(isinstance(rna_input_seq, dict))
        options['rna_input_seq'] = rna_input_seq

        smiles_input = options.get('smiles_input', {})
        assert(isinstance(smiles_input, dict))
        options['smiles_input'] = smiles_input

        check_reaction = options.get('check_reaction', True)
        assert(isinstance(check_reaction, bool))
        options['check_reaction'] = check_reaction

        gen_metabolites = options.get('gen_metabolites', True)
        assert(isinstance(gen_metabolites, bool))
        options['gen_metabolites'] = gen_metabolites       

        gen_dna = options.get('gen_dna', True)
        assert(isinstance(gen_dna, bool))
        options['gen_dna'] = gen_dna

        gen_transcripts = options.get('gen_transcripts', True)
        assert(isinstance(gen_transcripts, bool))
        options['gen_transcripts'] = gen_transcripts

        gen_protein = options.get('gen_protein', True)
        assert(isinstance(gen_protein, bool))
        options['gen_protein'] = gen_protein

        gen_complexes = options.get('gen_complexes', True)
        assert(isinstance(gen_complexes, bool))
        options['gen_complexes'] = gen_complexes

        gen_distribution_init_concentrations = options.get('gen_distribution_init_concentrations', True)
        assert(isinstance(gen_distribution_init_concentrations, bool))
        options['gen_distribution_init_concentrations'] = gen_distribution_init_concentrations

        gen_observables = options.get('gen_observables', True)
        assert(isinstance(gen_observables, bool))
        options['gen_observables'] = gen_observables

        gen_kb_reactions = options.get('gen_kb_reactions', True)
        assert(isinstance(gen_kb_reactions, bool))
        options['gen_kb_reactions'] = gen_kb_reactions

        gen_dfba_objective = options.get('gen_dfba_objective', False)
        assert(isinstance(gen_dfba_objective, bool))
        options['gen_dfba_objective'] = gen_dfba_objective

        gen_kb_rate_laws = options.get('gen_kb_rate_laws', True)
        assert(isinstance(gen_kb_rate_laws, bool))
        options['gen_kb_rate_laws'] = gen_kb_rate_laws

        gen_environment = options.get('gen_environment', True)
        assert(isinstance(gen_environment, bool))
        options['gen_environment'] = gen_environment

    def gen_taxon(self):
        """ Generate taxon for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model
        
        ncbi_taxa = taxoniq.Taxon(kb.cell.taxon)
        taxon_name = ncbi_taxa.scientific_name
        taxon_rank = ncbi_taxa.rank.name
        model_taxon = wc_lang.core.Taxon(id=str(kb.cell.taxon), name=taxon_name, 
            model=model, rank=wc_lang.core.TaxonRank[taxon_rank])

    def gen_compartments(self):
        """ Generate compartments for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        if kb.cell.parameters.get_one(id='cell_volume'):
            mean_cell_volume = kb.cell.parameters.get_one(id='cell_volume').value
        else:
            raise ValueError('The cell object does not have the parameter cell_volume')        
        
        culture_volume = self.options['culture_volume']
        cell_density = self.options['cell_density']
        membrane_density = self.options['membrane_density']

        for comp in kb.cell.compartments:

            c = model.compartments.get_or_create(
                    id=comp.id, name=comp.name)

            c.init_density = model.parameters.create(
                id='density_' + c.id,
                units=unit_registry.parse_units('g l^-1'))

            if comp.id=='e':
                c.biological_type = wc_ontology['WC:extracellular_compartment']
                c.init_density.value = 1000.
                c.init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], 
                    mean=culture_volume, std=0)                
                
            elif '_m' in comp.id:
                c.physical_type = wc_ontology['WC:membrane_compartment']
                c.init_density.value = membrane_density
                organelle_fraction = kb.cell.compartments.get_one(id=comp.id[:comp.id.index('_')]).volumetric_fraction              
                c.init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], 
                    mean=4.836E-09*(mean_cell_volume*organelle_fraction)**(2/3), std=0)                

            else:
                c.init_density.value = cell_density
                organelle_fraction = kb.cell.compartments.get_one(id=comp.id).volumetric_fraction
                c.init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], 
                    mean=mean_cell_volume*organelle_fraction - 4.836E-09*(mean_cell_volume*organelle_fraction)**(2/3), std=0)

            volume = model.functions.create(id='volume_' + c.id, units=unit_registry.parse_units('l'))
            volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
            assert error is None, str(error)           

    def gen_parameters(self):
        """ Generate parameters for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        Avogadro = model.parameters.create(id='Avogadro',
                                        type = None,
                                        value = scipy.constants.Avogadro,
                                        units = unit_registry.parse_units('molecule mol^-1'))       
      
        # Create parameters from kb
        for param in kb.cell.parameters:
            model_param = model.parameters.create(
                            id=param.id,
                            name=param.name,                            
                            value=param.value,
                            units=param.units)
            if 'K_m' in param.id:
                model_param.type = wc_ontology['WC:K_m']
            elif 'k_cat' in param.id:
                model_param.type = wc_ontology['WC:k_cat']
                model_param.units = unit_registry.parse_units('molecule^-1 s^-1')
            else:
                model_param.type = None

            if param.references:
                for ref in param.references:
                    ref_model = model.references.get_or_create(
                        __type=wc_lang.Reference,
                        author=ref.authors,
                        title=ref.title,
                        publication=ref.journal,
                        volume=ref.volume,
                        issue=ref.issue,
                        pages=ref.pages,
                        year=ref.year,
                        comments=ref.comments, 
                        type=wc_ontology['WC:article'])
                    if not ref_model.id:
                        ref_model.id = 'ref_'+str(len(model.references))
                    model_param.references.append(ref_model)

            if param.identifiers:
                for identifier in param.identifiers:
                    identifier_model = wc_lang.Identifier(
                        namespace=identifier.namespace, id=identifier.id)
                    model_param.identifiers.append(identifier_model)

        # Standardize the units of doubling time
        if model.parameters.get_one(id='mean_doubling_time'):
            model_doubling_time = model.parameters.get_one(id='mean_doubling_time')
        else:
            raise ValueError('The cell object does not have the parameter mean_doubling_time')

        expr = unit_registry.parse_expression(str(model_doubling_time.units))
        scale = expr.to(unit_registry.parse_units('second'))
        conversion_factor = scale.magnitude
        model_doubling_time.value *= conversion_factor
        model_doubling_time.units = unit_registry.parse_units('s')

    def global_vars_from_input(self):
        """ Populate global variable if input transcript sequences are provided in the options  
        """
        rna_input_seq = self.options['rna_input_seq']
        selenoproteome =self.options['selenoproteome']

        for Id, seq in rna_input_seq.items():
            gvar.transcript_ntp_usage[Id] = {
                'A': seq.upper().count('A'),
                'C': seq.upper().count('C'),
                'G': seq.upper().count('G'),
                'U': seq.upper().count('U'),
                'len': len(seq)
                }

    def gen_metabolites(self):
        """ Generate metabolites for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.core.MetaboliteSpeciesType)

        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)                    

    def gen_dna(self):
        """ Generate DNAs for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.core.DnaSpeciesType)

        for kb_species_type in kb_species_types:
            if 'M' in kb_species_type.id:
                self.gen_species_type(kb_species_type, ['m'])
            else:
                self.gen_species_type(kb_species_type, ['n'])

    def gen_transcripts(self):
        """ Generate transcripts (mature RNAs) for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.eukaryote.TranscriptSpeciesType)

        count = 0
        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)
            count += 1
            if count % 100 == 0:
                print('{}/{} of the transcripts have been initialized'.format(
                    count, len(kb_species_types)))

    def gen_protein(self):
        """ Generate proteins for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.eukaryote.ProteinSpeciesType)

        count = 0
        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)
            count += 1
            if count % 100 == 0:
                print('{}/{} of the proteins have been initialized'.format(
                    count, len(kb_species_types)))

    def gen_complexes(self):
        """ Generate complexes for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.core.ComplexSpeciesType)

        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)

    def gen_species_type(self, kb_species_type, extra_compartment_ids=None):
        """ Generate a model species type and species

        Args:
            kb_species_type (:obj:`wc_kb.SpeciesType`): knowledge base species type
            extra_compartment_ids (:obj:`list` of :obj:`str`, optional): compartment ids of
                additional species that should be created beyond those represented in the KB

        Returns:
            * :obj:`wc_lang.SpeciesType`: model species type
        """

        model = self.model
        ph = self.options['ph']
        selenoproteome = self.options['selenoproteome']        

        model_species_type = model.species_types.get_or_create(id=kb_species_type.id)
        model_species_type.name = kb_species_type.name
        model_species_type.structure = wc_lang.ChemicalStructure()
        model_species_type.comments = kb_species_type.comments

        if isinstance(kb_species_type, wc_kb.core.MetaboliteSpeciesType):
            model_species_type.type = wc_ontology['WC:metabolite'] # metabolite
            if kb_species_type.name == 'electron':
                model_species_type.structure.value = '[*-]'
                model_species_type.structure.format = wc_lang.ChemicalStructureFormat.SMILES
                model_species_type.structure.empirical_formula = EmpiricalFormula()
                model_species_type.structure.molecular_weight = 0.
                model_species_type.structure.charge = -1
            elif kb_species_type.name == 'proton' or kb_species_type.name == 'proton motive force':
                model_species_type.structure.value = '[1H+]'
                model_species_type.structure.format = wc_lang.ChemicalStructureFormat.SMILES
                model_species_type.structure.empirical_formula = EmpiricalFormula('H')
                model_species_type.structure.molecular_weight = 1.008
                model_species_type.structure.charge = 1
            elif kb_species_type.name == 'dihydrogen':
                model_species_type.structure.value = '[H][H]'
                model_species_type.structure.format = wc_lang.ChemicalStructureFormat.SMILES
                model_species_type.structure.empirical_formula = EmpiricalFormula('H2')
                model_species_type.structure.molecular_weight = 2.016
                model_species_type.structure.charge = 0    
            else:    
                inchi_str = kb_species_type.properties.get_one(property='structure')
                if inchi_str:
                    smiles, formula, charge, mol_wt = self.structure_to_smiles_and_props(
                        kb_species_type.id, inchi_str.get_value(), ph)
                    model_species_type.structure.value = smiles
                    model_species_type.structure.format = wc_lang.ChemicalStructureFormat.SMILES
                else:
                    formula = kb_species_type.get_empirical_formula()
                    charge = kb_species_type.get_charge()
                    mol_wt = kb_species_type.get_mol_wt()
                model_species_type.structure.empirical_formula = formula
                model_species_type.structure.molecular_weight = mol_wt
                model_species_type.structure.charge = charge 

        elif isinstance(kb_species_type, wc_kb.core.DnaSpeciesType):
            model_species_type.type = wc_ontology['WC:DNA'] # DNA
            model_species_type.structure.empirical_formula = kb_species_type.get_empirical_formula()
            model_species_type.structure.molecular_weight = kb_species_type.get_mol_wt()
            model_species_type.structure.charge = kb_species_type.get_charge()

        elif isinstance(kb_species_type, wc_kb.eukaryote.TranscriptSpeciesType):
            model_species_type.type = wc_ontology['WC:RNA'] # RNA
            if model_species_type.id in self.options['rna_input_seq']:
                seq = self.options['rna_input_seq'][model_species_type.id]
            else:
                seq = kb_species_type.get_seq()
                gvar.transcript_ntp_usage[model_species_type.id] = {
                    'A': seq.upper().count('A'),
                    'C': seq.upper().count('C'),
                    'G': seq.upper().count('G'),
                    'U': seq.upper().count('U'),
                    'len': len(seq)
                    }
            model_species_type.structure.empirical_formula = kb_species_type.get_empirical_formula(
                seq_input=seq)
            model_species_type.structure.molecular_weight = kb_species_type.get_mol_wt(
                seq_input=seq)
            model_species_type.structure.charge = kb_species_type.get_charge(
                seq_input=seq)

        elif isinstance(kb_species_type, wc_kb.eukaryote.ProteinSpeciesType):
            model_species_type.type = wc_ontology['WC:protein'] # protein
            table = 2 if 'M' in kb_species_type.transcript.gene.polymer.id else 1
            cds = self.options['cds']
            _, raw_seq, start_codon = kb_species_type.get_seq_and_start_codon(table=table, cds=cds)
            if kb_species_type.transcript.gene.id in selenoproteome:
                processed_seq = raw_seq[:-1] if raw_seq.endswith('*') else raw_seq
                protein_seq = ''.join(i if i!='*' else 'U' for i in processed_seq)
            else:                                            
                protein_seq = ''.join(i for i in raw_seq if i!='*')            
            self.populate_protein_aa_usage(model_species_type.id, protein_seq)
            gvar.protein_aa_usage[model_species_type.id]['start_aa'] = protein_seq[0]
            gvar.protein_aa_usage[model_species_type.id]['start_codon'] = str(start_codon).upper()                        
            _, _, _, determined = self.determine_protein_structure_from_aa(
                model_species_type.id, gvar.protein_aa_usage[model_species_type.id])
            if not determined:    
                model_species_type.structure.empirical_formula = kb_species_type.get_empirical_formula(
                    seq_input=protein_seq)
                model_species_type.structure.molecular_weight = kb_species_type.get_mol_wt(
                    seq_input=protein_seq)
                model_species_type.structure.charge = kb_species_type.get_charge(
                    seq_input=protein_seq)       

        elif isinstance(kb_species_type, wc_kb.core.ComplexSpeciesType):
            model_species_type.type = wc_ontology['WC:pseudo_species'] # pseudo specie
            formula = EmpiricalFormula()
            charge = 0
            weight = 0
            for subunit in kb_species_type.subunits:
                subunit_id = subunit.species_type.id
                coef = abs(max(1, subunit.coefficient))
                if isinstance(subunit.species_type, wc_kb.eukaryote.ProteinSpeciesType):
                    subunit_model_species_type = model.species_types.get_one(id=subunit_id)
                    if subunit_model_species_type:
                        formula += subunit_model_species_type.structure.empirical_formula * coef 
                        charge += subunit_model_species_type.structure.charge * coef
                        weight += subunit_model_species_type.structure.molecular_weight * coef
                    else:    
                        table = 2 if 'M' in subunit.species_type.transcript.gene.polymer.id else 1
                        cds = self.options['cds']
                        _, raw_seq, start_codon = subunit.species_type.get_seq_and_start_codon(table=table, cds=cds)
                        if subunit.species_type.transcript.gene.id in selenoproteome:
                            processed_seq = raw_seq[:-1] if raw_seq.endswith('*') else raw_seq
                            protein_seq = ''.join(i if i!='*' else 'U' for i in processed_seq)
                        else:                                            
                            protein_seq = ''.join(i for i in raw_seq if i!='*')
                        self.populate_protein_aa_usage(subunit_id, protein_seq)
                        gvar.protein_aa_usage[subunit_id]['start_aa'] = protein_seq[0]
                        gvar.protein_aa_usage[subunit_id]['start_codon'] = str(start_codon).upper()                        
                        sub_formula, sub_mol_wt, sub_charge, determined = \
                            self.determine_protein_structure_from_aa(
                            subunit_id, gvar.protein_aa_usage[subunit_id])
                        if not determined:
                            formula += subunit.species_type.get_empirical_formula(seq_input=protein_seq) * coef
                            charge += subunit.species_type.get_charge(seq_input=protein_seq) * coef
                            weight += subunit.species_type.get_mol_wt(seq_input=protein_seq) * coef
                        else:
                            formula += sub_formula * coef 
                            charge += sub_charge * coef
                            weight += sub_mol_wt * coef
                elif isinstance(subunit.species_type, wc_kb.eukaryote.TranscriptSpeciesType):
                    subunit_model_species_type = model.species_types.get_one(id=subunit_id)
                    if subunit_model_species_type:
                        formula += subunit_model_species_type.structure.empirical_formula * coef 
                        charge += subunit_model_species_type.structure.charge * coef
                        weight += subunit_model_species_type.structure.molecular_weight * coef
                    else:
                        if subunit_id not in gvar.transcript_ntp_usage:
                            seq = subunit.species_type.get_seq()
                            gvar.transcript_ntp_usage[subunit.species_type.id] = {
                                'A': seq.upper().count('A'),
                                'C': seq.upper().count('C'),
                                'G': seq.upper().count('G'),
                                'U': seq.upper().count('U'),
                                'len': len(seq)
                                }
                        else:
                            seq = self.options['rna_input_seq'][subunit_id]

                        formula += subunit.species_type.get_empirical_formula(seq_input=seq) * coef
                        weight += subunit.species_type.get_mol_wt(seq_input=seq) * coef
                        charge += subunit.species_type.get_charge(seq_input=seq) * coef

                else:
                    inchi_str = subunit.species_type.properties.get_one(property='structure')
                    if inchi_str:
                        _, sub_formula, sub_charge, sub_mol_wt = self.structure_to_smiles_and_props(
                            subunit.species_type.id, inchi_str.get_value(), ph)
                    else:
                        sub_formula = subunit.species_type.get_empirical_formula()
                        sub_charge = subunit.species_type.get_charge()
                        sub_mol_wt = subunit.species_type.get_mol_wt()    
                    
                    formula += sub_formula * coef
                    charge += sub_charge * coef
                    weight += sub_mol_wt * coef        
            
            model_species_type.structure.empirical_formula = formula
            model_species_type.structure.molecular_weight = weight
            model_species_type.structure.charge = charge

        else:
            raise ValueError('Unsupported species type: {}'.format(
                kb_species_type.__class__.__name__))

        # Create species
        if isinstance(kb_species_type, wc_kb.core.ComplexSpeciesType):
            subunit_compartments = [[s.compartment.id for s in sub.species_type.species]
                for sub in kb_species_type.subunits]
            if len(subunit_compartments) == 1:
                shared_compartments = set(subunit_compartments[0])
            else:
                shared_compartments = set([])
                for i in range(len(subunit_compartments)):
                    shared_compartments = (set(subunit_compartments[i])
                        if i==0 else shared_compartments).intersection(
                        set(subunit_compartments[i+1]) if i<(len(subunit_compartments)-1) else shared_compartments)
            # Combine compartments where all the subunits exist, where catalyzed reactions occur and the additionally defined extra
            compartment_ids = set(list(shared_compartments) + [s.compartment.id for s in kb_species_type.species] +
                              (extra_compartment_ids or []))
        else:
            compartment_ids = set([s.compartment.id for s in kb_species_type.species] +
                              (extra_compartment_ids or []))

        for compartment_id in compartment_ids:
            model_compartment = model.compartments.get_one(id=compartment_id)
            model_species = model.species.get_or_create(species_type=model_species_type, compartment=model_compartment)
            model_species.id = model_species.gen_id()

        return model_species_type

    def gen_distribution_init_concentrations(self):
        """ Generate concentrations for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model
        media = self.options['media']

        Avogadro = model.parameters.get_one(id='Avogadro')

        for conc in kb.cell.concentrations:
            species_comp_model = model.compartments.get_one(id=conc.species.compartment.id)

            species_type = model.species_types.get_or_create(id=conc.species.species_type.id)
            species = model.species.get_or_create(species_type=species_type, compartment=species_comp_model)
            species.id = species.gen_id()

            if conc.units == unit_registry.parse_units('molecule'):
                mean_concentration = conc.value
            elif conc.units == unit_registry.parse_units('M'):
                mean_concentration = conc.value * Avogadro.value * species_comp_model.init_volume.mean
            else:
                raise Exception('Unsupported units: {}'.format(conc.units.name))

            conc_model = model.distribution_init_concentrations.create(
                species=species,
                mean=mean_concentration,
                units=unit_registry.parse_units('molecule'),
                comments=conc.comments)
            conc_model.id = conc_model.gen_id()

            if conc.references:
                for ref in conc.references:
                    ref_model = model.references.get_or_create(
                        __type=wc_lang.Reference,
                        author=ref.authors,
                        title=ref.title,
                        publication=ref.journal,
                        volume=ref.volume,
                        issue=ref.issue,
                        pages=ref.pages,
                        year=ref.year,
                        comments=ref.comments, 
                        type=wc_ontology['WC:article'])
                    if not ref_model.id:
                        ref_model.id = 'ref_'+str(len(model.references))
                    conc_model.references.append(ref_model)

            if conc.identifiers:
                for identifier in conc.identifiers:
                    identifier_model = wc_lang.Identifier(
                        namespace=identifier.namespace, id=identifier.id)
                    conc_model.identifiers.append(identifier_model)

        for chromosome in kb.cell.species_types.get(__type=wc_kb.core.DnaSpeciesType):
            model_species_type = model.species_types.get_or_create(id=chromosome.id)
            model_species = model.species.get_or_create(species_type=model_species_type)
            conc_model = model.distribution_init_concentrations.create(
                species=model_species,
                mean=chromosome.ploidy, 
                units=unit_registry.parse_units('molecule'),
                )
            conc_model.id = conc_model.gen_id()

        for Id, (conc, refs, comments) in media.items():
            species_type = model.species_types.get_or_create(id=Id)
            species_comp_model = model.compartments.get_one(id='e')
            species = model.species.get_or_create(species_type=species_type, compartment=species_comp_model)
            species.id = species.gen_id()

            conc_model = model.distribution_init_concentrations.create(
                species=species,
                mean=conc * Avogadro.value * species_comp_model.init_volume.mean,
                units=unit_registry.parse_units('molecule'),
                comments=comments,
                )
            conc_model.id = conc_model.gen_id()

            if refs:
                for ref in refs:
                    ref_model = model.references.get_one(title=ref.title)
                    if ref_model:
                        conc_model.references.append(ref_model)
                    else:    
                        ref.model = model
                        ref.id = 'ref_'+str(len(model.references))
                        conc_model.references.append(ref)
            
    def gen_observables(self):
        """ Generate observables for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        for kb_observable in kb.cell.observables:
            all_species = {}
            all_observables = {}

            for kb_species in kb_observable.expression.species:
                kb_species_type = kb_species.species_type
                kb_compartment = kb_species.compartment
                model_species_type = model.species_types.get_one(
                    id=kb_species_type.id)
                model_species = model_species_type.species.get_one(
                    compartment=model.compartments.get_one(id=kb_compartment.id))
                all_species[model_species.gen_id()] = model_species

            for kb_observable_observable in kb_observable.expression.observables:
                model_observable_observable = model.observables.get_or_create(
                    id=kb_observable_observable.id)
                all_observables[model_observable_observable.id] = model_observable_observable

            model_observable_expression, error = wc_lang.ObservableExpression.deserialize(
                kb_observable.expression.expression, {
                wc_lang.Species: all_species,
                wc_lang.Observable: all_observables,
                })
            assert error is None, str(error)

            model_observable = model.observables.get_or_create(
                id=kb_observable.id,
                name=kb_observable.name,
                expression=model_observable_expression)

    def gen_kb_reactions(self):
        """ Generate reactions encoded within KB
        """

        kb = self.knowledge_base
        model = self.model
        submodel = model.submodels.get_or_create(id='metabolism', 
            framework=wc_ontology['WC:dynamic_flux_balance_analysis'])              
 
        for kb_rxn in kb.cell.reactions:
            if len(kb_rxn.participants)==1 and kb_rxn.participants[0].species.compartment.id=='e' and \
                not model.distribution_init_concentrations.get_one(
                species=model.species.get_one(id=kb_rxn.participants[0].species.id())):
                pass
            else:    
                model_rxn = model.reactions.create(
                    submodel=submodel,
                    id=kb_rxn.id + '_kb',
                    name=kb_rxn.name,
                    reversible=kb_rxn.reversible,
                    comments=kb_rxn.comments)

                delta_formula = EmpiricalFormula()
                delta_charge = 0.
                proton_participant = 0
                for participant in kb_rxn.participants:
                    kb_species = participant.species
                    model_species_type = model.species_types.get_one(id=kb_species.species_type.id)
                    model_compartment = model.compartments.get_one(id=kb_species.compartment.id)
                    model_species = model_species_type.species.get_one(compartment=model_compartment)
                    
                    if model_species_type.name=='proton':
                        proton_participant += 1

                    if model_species is None:
                        model_species = self.gen_species_type(kb_species.species_type, model_compartment)

                    model_rxn.participants.add(
                        model_species.species_coefficients.get_or_create(coefficient=participant.coefficient))

                    # Check element and charge balance          
                    delta_formula += model_species_type.structure.empirical_formula * participant.coefficient
                    delta_charge += model_species_type.structure.charge * participant.coefficient

                # Correct proton and charge balance at the pH at which metabolite properties are determined
                if self.options['check_reaction']:
                    if delta_charge and len(delta_formula)==1 and delta_charge==delta_formula['H']:
                        proton_species_type = model.species_types.get_one(name='proton')
                        if proton_participant:
                            proton_coef = [i for i in model_rxn.participants if i.species.species_type==proton_species_type][0]
                            model_rxn.participants.discard(proton_coef)
                            model_rxn.participants.add(
                                proton_coef.species.species_coefficients.get_or_create(coefficient=proton_coef.coefficient-delta_charge))
                        else:
                            comp_id = set([i.species.compartment.id for i in model_rxn.participants]).pop()
                            proton_compartment = model.compartments.get_one(id=comp_id)
                            proton_species = model.species.get_or_create(species_type=proton_species_type, compartment=proton_compartment)
                            proton_species.id = proton_species.gen_id()
                            model_rxn.participants.add(
                                proton_species.species_coefficients.get_or_create(coefficient=-delta_charge))

        # For testing purpose
        if self.options['gen_dfba_objective']:
            submodel.dfba_obj = wc_lang.DfbaObjective(model=model)
            submodel.dfba_obj.id = submodel.dfba_obj.gen_id()
            obj_expression = model_rxn.id
            dfba_obj_expression, error = wc_lang.DfbaObjectiveExpression.deserialize(
                obj_expression, {wc_lang.Reaction: {model_rxn.id: model_rxn}})
            assert error is None, str(error)
            submodel.dfba_obj.expression = dfba_obj_expression            

    def gen_kb_rate_laws(self):
        """ Generate the rate laws for reactions encoded in the knowledge base """
        kb = self.knowledge_base
        model = self.model

        Avogadro = model.parameters.get_or_create(id='Avogadro')

        for kb_rxn in kb.cell.reactions:

            model_rxn = model.reactions.get_one(id=kb_rxn.id + '_kb')

            for kb_rate_law in kb_rxn.rate_laws:
                all_parameters = {}
                all_parameters[Avogadro.id] = Avogadro
                all_species = {}
                all_observables = {}
                all_volumes = {}
                
                model_expression = kb_rate_law.expression.expression

                for observable in kb_rate_law.expression.observables:
                    all_observables[observable.id] = model.observables.get_one(id=observable.id)

                for species in kb_rate_law.expression.species:
                    model_species_type = model.species_types.get_one(id=species.species_type.id)
                    model_compartment = model.compartments.get_one(id=species.compartment.id)
                    if model_compartment.init_density.function_expressions:
                        volume = model_compartment.init_density.function_expressions[0].function
                        all_volumes[volume.id] = volume
                    model_species = model_species_type.species.get_one(compartment=model_compartment)
                    all_species[model_species.gen_id()] = model_species                    

                for param in kb_rate_law.expression.parameters:
                    all_parameters[param.id] = model.parameters.get_one(id=param.id)
                    if 'K_m' in param.id:
                        model_compartment = model.compartments.get_one(id=param.id.split('_')[-1])
                        if model_compartment.init_density.function_expressions:
                            volume = model_compartment.init_density.function_expressions[0].function
                            unit_adjusted_term = '{} * {} * {}'.format(param.id, Avogadro.id, volume.id)
                        else:
                            volume = model_compartment.init_volume.mean    
                            unit_adjusted_term = '{} * {} * {}'.format(param.id, Avogadro.id, volume)
                        model_expression = model_expression.replace(param.id, unit_adjusted_term)            

                rate_law_expression, error = wc_lang.RateLawExpression.deserialize(
                    model_expression, {
                    wc_lang.Parameter: all_parameters,
                    wc_lang.Species: all_species,
                    wc_lang.Observable: all_observables,
                    wc_lang.Function: all_volumes,
                    })
                assert error is None, str(error)

                model_rate_law = model.rate_laws.create(
                    expression=rate_law_expression,
                    reaction=model_rxn,
                    direction=wc_lang.RateLawDirection[kb_rate_law.direction.name],
                    comments=kb_rate_law.comments)
                model_rate_law.id = model_rate_law.gen_id()

    def gen_environment(self):            
        """ Generate the environment, i.e. temperature, for the simulated cells """

        model = self.model
        environment = self.options['environment']

        if environment:
            wc_lang.Environment(
                id=environment['id'],
                name=environment['name'],
                model=model,
                temp=environment['temperature'],
                comments=environment['comments'])

    def structure_to_smiles_and_props(self, met_id, structure, ph):
        """ Convert InChI or SMILES string in the knowledge base 
            to a SMILES string at specific pH and calculate properties such
            as empirical formula, charge and molecular weight 

        Args:
            met_id (:obj:`str`): id of metabolite
            structure (:obj:`str`): InChI or SMILES string
            ph (:obj:`float`): pH at which the properties should be determined

        Returns:
            :obj:`str`: SMILES string
            :obj:`wc_utils.util.chem.core.EmpiricalFormula`: empirical formula
            :obj:`int`: charge
            :obj:`float`: molecular weight    
        """
        smiles_input = self.options['smiles_input']
        
        if met_id in smiles_input:
            smiles = smiles_input[met_id]
        else:    
            structure_type = 'inchi' if 'InChI=' in structure else 'smiles'
            smiles = get_major_micro_species(structure, structure_type, 'smiles', ph=ph)        
        
        mol = openbabel.OBMol()
        conv = openbabel.OBConversion()
        conv.SetInFormat('smi')
        conv.SetOptions('c', conv.OUTOPTIONS)
        conv.ReadString(mol, smiles)        
        empirical_formula = OpenBabelUtils.get_formula(mol)
        charge = mol.GetTotalCharge()
        mol_wt = empirical_formula.get_molecular_weight()
        
        return smiles, empirical_formula, charge, mol_wt

    def populate_protein_aa_usage(self, protein_id, seq):
        """ Populate a global variable dictionary of amino acid
            usage in a protein given its sequence

        Args:
            protein_id (:obj:`str`): protein ID
            seq (:obj:`Bio.Seq.Seq`): sequence    
        """
        gvar.protein_aa_usage[protein_id] = {
                'len': len(seq) - seq.count('*'),
                '*': seq.count('*'),  # Symbol used in Bio.Seq.Seq when cds is set to False  
                'A': seq.count('A'),  # Ala: Alanine (C3 H7 N O2)
                'R': seq.count('R'),  # Arg: Arginine (C6 H14 N4 O2)
                'N': seq.count('N'),  # Asn: Asparagine (C4 H8 N2 O3)
                'D': seq.count('D'),  # Asp: Aspartic acid (C4 H7 N O4)
                'C': seq.count('C'),  # Cys: Cysteine (C3 H7 N O2 S)
                'Q': seq.count('Q'),  # Gln: Glutamine (C5 H10 N2 O3)
                'E': seq.count('E'),  # Glu: Glutamic acid (C5 H9 N O4)
                'G': seq.count('G'),  # Gly: Glycine (C2 H5 N O2)
                'H': seq.count('H'),  # His: Histidine (C6 H9 N3 O2)
                'I': seq.count('I'),  # Ile: Isoleucine (C6 H13 N O2)
                'L': seq.count('L'),  # Leu: Leucine (C6 H13 N O2)
                'K': seq.count('K'),  # Lys: Lysine (C6 H14 N2 O2)
                'M': seq.count('M'),  # Met: Methionine (C5 H11 N O2 S)
                'F': seq.count('F'),  # Phe: Phenylalanine (C9 H11 N O2)
                'P': seq.count('P'),  # Pro: Proline (C5 H9 N O2)
                'S': seq.count('S'),  # Ser: Serine (C3 H7 N O3)
                'T': seq.count('T'),  # Thr: Threonine (C4 H9 N O3)
                'W': seq.count('W'),  # Trp: Tryptophan (C11 H12 N2 O2)
                'Y': seq.count('Y'),  # Tyr: Tyrosine (C9 H11 N O3)
                'V': seq.count('V'),  # Val: Valine (C5 H11 N O2)
                'U': seq.count('U'),  # Selcys: Selenocysteine (C3 H7 N O2 Se)
                'start_aa': seq[0] if seq[0]!='*' else seq[1] # Amino acid at start codon
            }

    def determine_protein_structure_from_aa(self, 
            polymer_id, count):
        """ Determine the empirical formula, molecular weight and charge of
            a protein based on the structural information of its metabolite
            amino acid monomers to ensure consistency with the pH

        Args:
            polymer_id (:obj:`str`): polymer ID
            count (:obj:`dict`): dictionary showing the count of each amino
                acid in the protein

        Returns:
            :obj:`wc_utils.util.chem.EmpiricalFormula`: protein empirical formula
            :obj:`float`: protein molecular weight
            :obj:`int`: protein charge
            :obj:`bool`: True if protein structure has been successfully determined
                from the metabolite monomer, else False             
        """                
        model = self.model
        amino_acid_id_conversion = self.options['amino_acid_id_conversion']
        
        total_empirical_formula = EmpiricalFormula()
        total_molecular_weight = 0
        total_charge = 0
        polymer_len = 0

        if not amino_acid_id_conversion:            
            return EmpiricalFormula(), 0, 0, False        
        
        for standard_id, met_id in amino_acid_id_conversion.items():
            if standard_id in count:
                monomer_species_type = model.species_types.get_one(id=met_id)
                if monomer_species_type:
                    total_empirical_formula += \
                        monomer_species_type.structure.empirical_formula * count[standard_id]
                    total_molecular_weight += \
                        monomer_species_type.structure.molecular_weight * count[standard_id]
                    total_charge += monomer_species_type.structure.charge * count[standard_id]
                    polymer_len += count[standard_id]    
                else:
                    return EmpiricalFormula(), 0, 0, False
            else:
                return EmpiricalFormula(), 0, 0, False        
        
        protein_empirical_formula = total_empirical_formula - \
            EmpiricalFormula('H{}O{}'.format(2 * (polymer_len - 1), polymer_len - 1))
        protein_molecular_weight = total_molecular_weight - \
            2 * (polymer_len - 1) * mendeleev.element('H').atomic_weight - \
            (polymer_len - 1) * mendeleev.element('O').atomic_weight 
        protein_charge = total_charge    
            
        model_species_type = model.species_types.get_one(id=polymer_id)
        if model_species_type:
            model_species_type.structure.empirical_formula = protein_empirical_formula
            model_species_type.structure.molecular_weight = protein_molecular_weight
            model_species_type.structure.charge = protein_charge
        
        return protein_empirical_formula, protein_molecular_weight, protein_charge, True        
