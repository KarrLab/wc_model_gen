""" Initalize the construction of wc_lang-encoded models from wc_kb-encoded knowledge base.

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-09
:Copyright: 2019, Karr Lab
:License: MIT
"""

from wc_utils.util.chem import EmpiricalFormula
from wc_utils.util.ontology import wcm_ontology
from wc_utils.util.units import unit_registry, are_units_equivalent
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

        self.gen_compartments()
        self.gen_parameters()

        if options['gen_dna']:
            self.gen_dna()

        if options['gen_pre_rnas']:
            self.gen_pre_rnas()

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

        cell_density = options.get('cell_density', 1040.)
        assert(isinstance(cell_density, float))
        options['cell_density'] = cell_density

        gen_dna = options.get('gen_dna', True)
        assert(isinstance(gen_dna, bool))
        options['gen_dna'] = gen_dna

        gen_pre_rnas = options.get('gen_pre_rnas', True)
        assert(isinstance(gen_pre_rnas, bool))
        options['gen_pre_rnas'] = gen_pre_rnas

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

    def gen_compartments(self):
        
        kb = self.knowledge_base
        model = self.model

        if kb.cell.properties.get_one(id='cell_volume'):
            mean_cell_volume = kb.cell.properties.get_one(id='cell_volume').value
        else:
            raise ValueError('The cell object does not have the property cell_volume')        
        
        cell_density = self.options['cell_density']

        for comp in kb.cell.compartments:
            c = model.compartments.get_or_create(
                    id=comp.id, 
                    name=comp.name, 
                    mean_init_volume=1E5*mean_cell_volume if comp.id=='e' else comp.volumetric_fraction*mean_cell_volume,
                    )
            c.init_density = model.parameters.create(
                id='density_' + c.id, 
                value=1000 if comp.id=='e' else cell_density, 
                units=unit_registry.parse_units('g l^-1'))
            volume = model.functions.create(id='volume_' + c.id, units=unit_registry.parse_units('l'))
            volume.expression, error = wc_lang.FunctionExpression.deserialize(f'{c.id} / {c.init_density.id}', {
                wc_lang.Compartment: {c.id: c},
                wc_lang.Parameter: {c.init_density.id: c.init_density},
                })
            assert error is None, str(error)

    def gen_parameters(self):
        kb = self.knowledge_base
        model = self.model

        # Create parameters out of properties
        if kb.cell.properties.get_one(id='mean_doubling_time'):
            doubling_time_kb = kb.cell.properties.get_one(id='mean_doubling_time')
        else:
            raise ValueError('The cell object does not have the property mean_doubling_time')

        if not isinstance(doubling_time_kb.units, unit_registry.Unit):
            ValueError('Unsupported units "{}"'.format(doubling_time_kb.units))

        expr = unit_registry.parse_expression(str(doubling_time_kb.units))
        scale = expr.to(unit_registry.parse_units('second'))
        conversion_factor = scale.magnitude

        model.parameters.get_or_create(id='mean_doubling_time',
                                       type=None,
                                       value=doubling_time_kb.value * conversion_factor,
                                       units=unit_registry.parse_units('s'))
 
        # Create parameters from kb
        for param in kb.cell.parameters:
            model_param = model.parameters.create(
                            id=param.id,                            
                            value=param.value,
                            units=param.units)
            if 'K_m' in param.id:
                model_param.type = wcm_ontology['WCM:K_m']
            elif 'k_cat' in param.id:
                model_param.type = wcm_ontology['WCM:k_cat']
            else:
                model_param.type = None

    def gen_dna(self):
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.core.DnaSpeciesType)

        for kb_species_type in kb_species_types:
            if 'MT' in kb_species_type.id or 'mt' in kb_species_type.id:
                self.gen_species_type(kb_species_type, ['m'])
            else:    
                self.gen_species_type(kb_species_type, ['n'])                

    def gen_pre_rnas(self):
        """ Generate pre-RNAs for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.eukaryote_schema.PreRnaSpeciesType)

        for kb_species_type in kb_species_types:
            chromosome = kb_species_type.gene.polymer
            if 'MT' in chromosome.id or 'mt' in chromosome.id:
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
            model_species_type.type = wcm_ontology['WCM:DNA'] # DNA

        elif isinstance(kb_species_type, wc_kb.eukaryote_schema.PreRnaSpeciesType):
            model_species_type.type = wcm_ontology['WCM:pseudo_species'] # pseudo species

        elif isinstance(kb_species_type, wc_kb.eukaryote_schema.TranscriptSpeciesType):
            model_species_type.type = wcm_ontology['WCM:RNA'] # RNA

        elif isinstance(kb_species_type, wc_kb.eukaryote_schema.ProteinSpeciesType):
            model_species_type.type = wcm_ontology['WCM:protein'] # protein

        elif isinstance(kb_species_type, wc_kb.core.MetaboliteSpeciesType):
            model_species_type.type = wcm_ontology['WCM:metabolite'] # metabolite
            model_species_type.structure = kb_species_type.structure    
            
        elif isinstance(kb_species_type, wc_kb.core.ComplexSpeciesType):
            model_species_type.type = wcm_ontology['WCM:pseudo_species'] # pseudo species

        else:
            raise ValueError('Unsupported species type: {}'.format(
                kb_species_type.__class__.__name__))

        if kb_species_type.get_empirical_formula():
            model_species_type.empirical_formula = EmpiricalFormula(kb_species_type.get_empirical_formula())

        model_species_type.molecular_weight = kb_species_type.get_mol_wt()
        model_species_type.charge = kb_species_type.get_charge()
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

        for conc in kb.cell.concentrations:
            species_comp_model = model.compartments.get_one(id=conc.species.compartment.id)

            species_type = model.species_types.get_or_create(id=conc.species.species_type.id)
            species = model.species.get_or_create(species_type=species_type, compartment=species_comp_model)
            species.id = species.gen_id()

            conc_model = model.distribution_init_concentrations.create(
                species=species,
                mean=conc.value, units=unit_registry.parse_units('M'),
                comments=conc.comments)
            conc_model.id = conc_model.gen_id()

            if conc.references:
                for ref in conc.references:
                    ref_model = wc_lang.Reference(model=model, id=ref.id, name=ref.standard_id, 
                        type=wcm_ontology['WCM:article'])
                    conc_model.references.append(ref_model)

            if conc.database_references:
                for db_ref in conc.database_references:
                    db_ref_model = wc_lang.DatabaseReference(model=model, 
                        database=db_ref.database, id=db_ref.id)
                    conc_model.db_refs.append(db_ref_model)                

    def gen_observables(self):
        """ Generate observables for the model from knowledge base """
        kb = self.knowledge_base
        model = self.model
        
        for kb_observable in kb.cell.observables:
            model_observable = model.observables.get_or_create(id=kb_observable.id)            
            model_observable.name = kb_observable.name
            model_observable.expression = wc_lang.ObservableExpression(
                expression=kb_observable.expression.expression)
            
            for kb_species in kb_observable.expression.species:
                kb_species_type = kb_species.species_type
                kb_compartment = kb_species.compartment
                model_species_type = model.species_types.get_one(
                    id=kb_species_type.id)
                model_species = model_species_type.species.get_one(
                    compartment=model.compartments.get_one(id=kb_compartment.id))
                model_observable.expression.species.append(model_species)
                
            for kb_observable_observable in kb_observable.expression.observables:
                model_observable_observable = model.observables.get_or_create(
                    id=kb_observable_observable.id)
                model_observable.expression.observables.append(model_observable_observable)

    def gen_kb_reactions(self):
        """ Generate the reactions encoded within the knowledge base """
        kb = self.knowledge_base
        model = self.model

        for kb_rxn in kb.cell.reactions:
            submodel_id = kb_rxn.submodel
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
        Avogadro.type = None
        Avogadro.value = scipy.constants.Avogadro
        Avogadro.units = unit_registry.parse_units('molecule mol^-1')

        for kb_rxn in kb.cell.reactions:

            model_rxn = model.reactions.get_one(id=kb_rxn.id + '_kb')

            for kb_rate_law in kb_rxn.rate_laws:
                model_rate_law = model.rate_laws.create(
                    expression= wc_lang.RateLawExpression(expression=kb_rate_law.expression.expression),
                    reaction=model_rxn,
                    direction=wc_lang.RateLawDirection[kb_rate_law.direction.name],
                    comments=kb_rate_law.comments)
                model_rate_law.id = model_rate_law.gen_id()

                for param in kb_rate_law.expression.parameters:
                    model_rate_law.expression.parameters.append(
                        model.parameters.get_one(id=param.id))
                    if 'K_m' in param.id:
                        volume = model.compartments.get_one(id=param.id[param.id.rfind('_')+1:]).init_density.function_expressions[0].function
                        unit_adjusted_term = '{} * {} * {}'.format(param.id, Avogadro.id, volume.id)
                        model_rate_law.expression.expression = model_rate_law.expression.expression.replace(param.id, unit_adjusted_term)
                        
                for species in kb_rate_law.expression.species:
                    model_species_type = model.species_types.get_one(id=species.species_type.id)
                    model_compartment = model.compartments.get_one(id=species.compartment.id)                    
                    model_rate_law.expression.species.append(
                        model_species_type.species.get_one(compartment=model_compartment))

                for observable in kb_rate_law.expression.observables:
                    model_rate_law.expression.observables.append(
                        model.observables.get_one(id=observable.id))
