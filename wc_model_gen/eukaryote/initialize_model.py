""" Initialize the construction of wc_lang-encoded models from wc_kb-encoded knowledge base.

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-09
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.chem import EmpiricalFormula
from wc_onto import onto as wc_ontology
from wc_utils.util.units import unit_registry
import ete3
import math
import numpy
import scipy.constants
import wc_kb
import wc_lang
import wc_model_gen


class InitializeModel(wc_model_gen.ModelComponentGenerator):
    """ Generate compartments """

    def run(self):
        self.clean_and_validate_options()
        options = self.options

        self.gen_taxon()
        self.gen_compartments()
        self.gen_parameters()

        if options['gen_dna']:
            self.gen_dna()

        if options['gen_transcripts']:
            self.gen_transcripts()

        if options['gen_protein']:
            self.gen_protein()

        if options['gen_metabolites']:
            self.gen_metabolites()

        if options['gen_complexes']:
            self.gen_complexes()

        if options['gen_distribution_init_concentrations']:
            self.gen_distribution_init_concentrations()

        if options['gen_observables']:
            self.gen_observables()

        if options['gen_kb_reactions']:
            self.gen_kb_reactions()

        if options['gen_kb_rate_laws']:
            self.gen_kb_rate_laws()

    def clean_and_validate_options(self):
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

        gen_dna = options.get('gen_dna', True)
        assert(isinstance(gen_dna, bool))
        options['gen_dna'] = gen_dna

        gen_transcripts = options.get('gen_transcripts', True)
        assert(isinstance(gen_transcripts, bool))
        options['gen_transcripts'] = gen_transcripts

        gen_protein = options.get('gen_protein', True)
        assert(isinstance(gen_protein, bool))
        options['gen_protein'] = gen_protein

        gen_metabolites = options.get('gen_metabolites', True)
        assert(isinstance(gen_metabolites, bool))
        options['gen_metabolites'] = gen_metabolites

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

        gen_kb_rate_laws = options.get('gen_kb_rate_laws', True)
        assert(isinstance(gen_kb_rate_laws, bool))
        options['gen_kb_rate_laws'] = gen_kb_rate_laws

    def gen_taxon(self):

        kb = self.knowledge_base
        model = self.model

        ncbi_taxa = ete3.NCBITaxa()
        taxon_name = ncbi_taxa.get_taxid_translator([kb.cell.taxon])[kb.cell.taxon]
        taxon_rank = ncbi_taxa.get_rank([kb.cell.taxon])[kb.cell.taxon]
        model_taxon = wc_lang.core.Taxon(id=str(kb.cell.taxon), name=taxon_name, model=model, 
            rank=wc_lang.core.TaxonRank[taxon_rank]) 

    def gen_compartments(self):

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
                c.init_density.value = 1000.
                c.init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], 
                    mean=culture_volume, std=0)                
                
            elif '_m' in comp.id:
                c.init_density.value = membrane_density
                organelle_fraction = kb.cell.compartments.get_one(id=comp.id[:comp.id.index('_')]).volumetric_fraction              
                c.init_volume = wc_lang.core.InitVolume(distribution=wc_ontology['WC:normal_distribution'], 
                    mean=4.836E-09*(mean_cell_volume*organelle_fraction)**(2/3), std=0)
                volume = model.functions.create(id='volume_' + c.id, units=unit_registry.parse_units('l'))                    
                volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                    wc_lang.Compartment: {c.id: c},
                    wc_lang.Parameter: {c.init_density.id: c.init_density},
                    })
                assert error is None, str(error)

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
                            value=param.value,
                            units=param.units)
            if 'K_m' in param.id:
                model_param.type = wc_ontology['WC:K_m']
            elif 'k_cat' in param.id:
                model_param.type = wc_ontology['WC:k_cat']
            else:
                model_param.type = None

            if param.references:
                for ref in param.references:
                    ref_model = wc_lang.Reference(model=model, id=ref.id, 
                        author=ref.authors,
                        title=ref.title,
                        publication=ref.journal,
                        volume=ref.volume,
                        issue=ref.issue,
                        pages=ref.pages,
                        year=ref.year,
                        comments=ref.comments, 
                        type=wc_ontology['WC:article'])
                    model_param.references.append(ref_model)

            if param.identifiers:
                for identifier in param.identifiers:
                    identifier_model = wc_lang.Identifier(model=model, 
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

    def gen_dna(self):
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
            __type=wc_kb.eukaryote_schema.TranscriptSpeciesType)

        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)

    def gen_protein(self):
        """ Generate proteins for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.eukaryote_schema.ProteinSpeciesType)

        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)

    def gen_metabolites(self):
        """ Generate metabolites for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.core.MetaboliteSpeciesType)

        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)

    def gen_complexes(self):
        """ Generate complexes for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.core.ComplexSpeciesType)

        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)

    def gen_species_type(self, kb_species_type, extra_compartment_ids=None):
        """ Generate a model species type

        Args:
            kb_species_type (:obj:`wc_kb.SpeciesType`): knowledge base species type
            extra_compartment_ids (:obj:`list` of :obj:`str`, optional): compartment ids of
                additional species that should be created beyond those represented in the KB

        Returns:
            * :obj:`wc_lang.SpeciesType`: model species type
        """

        model = self.model
        model_species_type = model.species_types.get_or_create(id=kb_species_type.id)
        model_species_type.name = kb_species_type.name

        if isinstance(kb_species_type, wc_kb.core.DnaSpeciesType):
            model_species_type.type = wc_ontology['WC:DNA'] # DNA

        elif isinstance(kb_species_type, wc_kb.eukaryote_schema.TranscriptSpeciesType):
            model_species_type.type = wc_ontology['WC:RNA'] # RNA

        elif isinstance(kb_species_type, wc_kb.eukaryote_schema.ProteinSpeciesType):
            model_species_type.type = wc_ontology['WC:protein'] # protein

        elif isinstance(kb_species_type, wc_kb.core.MetaboliteSpeciesType):
            model_species_type.type = wc_ontology['WC:metabolite'] # metabolite

        elif isinstance(kb_species_type, wc_kb.core.ComplexSpeciesType):
            model_species_type.type = wc_ontology['WC:pseudo_species'] # pseudo specie

        else:
            raise ValueError('Unsupported species type: {}'.format(
                kb_species_type.__class__.__name__))

        model_species_type.structure = wc_lang.ChemicalStructure()

        inchi_str = kb_species_type.properties.get_one(property='structure')
        if inchi_str:
            model_species_type.structure.value = inchi_str.get_value()

        formula = kb_species_type.get_empirical_formula()
        if formula:
            model_species_type.structure.empirical_formula = formula

        model_species_type.structure.molecular_weight = kb_species_type.get_mol_wt()
        model_species_type.structure.charge = kb_species_type.get_charge()
        model_species_type.comments = kb_species_type.comments

        if isinstance(kb_species_type, wc_kb.core.ComplexSpeciesType):
            subunit_compartments = [[s.compartment.id for s in sub.species_type.species]
                for sub in kb_species_type.subunits]

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
                    ref_model = wc_lang.Reference(model=model, id=ref.id, 
                        author=ref.authors,
                        title=ref.title,
                        publication=ref.journal,
                        volume=ref.volume,
                        issue=ref.issue,
                        pages=ref.pages,
                        year=ref.year,
                        comments=ref.comments, 
                        type=wc_ontology['WC:article'])
                    conc_model.references.append(ref_model)

            if conc.identifiers:
                for identifier in conc.identifiers:
                    identifier_model = wc_lang.Identifier(model=model, 
                        namespace=identifier.namespace, id=identifier.id)
                    conc_model.identifiers.append(identifier_model)

        for chromosome in kb.cell.species_types.get(__type=wc_kb.core.DnaSpeciesType):
            model_species_type = model.species_types.get_one(id=chromosome.id)
            model_species = model.species.get(species_type=model_species_type)[0]
            conc_model = model.distribution_init_concentrations.create(
                species=model_species,
                mean=chromosome.ploidy, 
                units=unit_registry.parse_units('molecule'),
                )
            conc_model.id = conc_model.gen_id()              

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
            TODO: rxn.submodel attribute is removed from KB, so submodel assignement needs to be taken care of here.
                  Since all rxns explicitly encoded in KB are metabolic ones, it is manually set to be metablism atm, but
                  not more sophisticated approach

        """

        kb = self.knowledge_base
        model = self.model

        for kb_rxn in kb.cell.reactions:
            submodel_id = 'Metabolism'
            submodel = model.submodels.get_or_create(id=submodel_id)

            model_rxn = model.reactions.create(
                submodel=submodel,
                id=kb_rxn.id + '_kb',
                name=kb_rxn.name,
                reversible=kb_rxn.reversible,
                comments=kb_rxn.comments)

            for participant in kb_rxn.participants:
                kb_species = participant.species
                model_species_type = model.species_types.get_one(id=kb_species.species_type.id)
                model_compartment = model.compartments.get_one(id=kb_species.compartment.id)
                model_species = model_species_type.species.get_one(compartment=model_compartment)

                if model_species is None:
                    model_species = self.gen_species_type(kb_species.species_type, model_compartment)

                model_rxn.participants.add(
                    model_species.species_coefficients.get_or_create(coefficient=participant.coefficient))

        # Generate complexation reactions
        for compl in kb.cell.species_types.get(__type=wc_kb.core.ComplexSpeciesType):
            model_compl_species_type = model.species_types.get_one(id=compl.id)

            for model_compl_species in model_compl_species_type.species:

                if all(model.species_types.get_one(id=subunit.species_type.id).species.get_one(
                    compartment=model_compl_species.compartment)!=None for subunit in compl.subunits):

                    submodel_id = 'Complexation'
                    submodel = model.submodels.get_or_create(id=submodel_id)

                    model_rxn = model.reactions.create(
                        submodel=submodel,
                        id=compl.id + '_' + model_compl_species.compartment.id,
                        name='Complexation of ' + compl.id + ' in ' + model_compl_species.compartment.name,
                        reversible=True)

                    for subunit in compl.subunits:
                        model_subunit_species = model.species_types.get_one(
                            id=subunit.species_type.id).species.get_one(compartment=model_compl_species.compartment)
                        model_rxn.participants.add(
                            model_subunit_species.species_coefficients.get_or_create(coefficient=-subunit.coefficient))

                    model_rxn.participants.add(
                            model_compl_species.species_coefficients.get_or_create(coefficient=1))

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
                        model_compartment = model.compartments.get_one(id=param.id[param.id.rfind('_')+1:])
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
