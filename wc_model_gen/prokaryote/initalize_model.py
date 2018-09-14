""" Initalize the construxction of wc_lang formatted models from wc_kb knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
- read cytosol volume from DB; currently there is onyl fractional volume?!
"""

from wc_lang import Species, Observable, ExpressionMethods
import wc_model_gen
import wc_lang
import wc_kb

class InitalizeModel(wc_model_gen.ModelComponentGenerator):
    """ Generate compartments """

    def run(self):
        self.gen_compartments()
        self.gen_parameters()
        self.clean_and_validate_options()
        options = self.options

        if options['gen_metabolic_species']:
            self.gen_metabolic_species()

        if options['gen_rna']:
            self.gen_rna()

        if options['gen_protein']:
            self.gen_protein()

        if options['gen_complexes']:
            self.gen_complexes()

        if options['gen_concentrations']:
            self.gen_concentrations()    

        if options['gen_observables']:
            self.gen_observables()

        if options['gen_kb_reactions']:
            self.gen_kb_reactions()

        if options['gen_kb_rate_laws']:
            self.gen_kb_rate_laws()

    def clean_and_validate_options(self):
        options = self.options

        gen_metabolic_species = options.get('gen_metabolic_species', True)
        assert(isinstance(gen_metabolic_species,bool))
        options['gen_metabolic_species'] = gen_metabolic_species

        gen_rna = options.get('gen_rna', True)
        assert(isinstance(gen_rna,bool))
        options['gen_rna'] = gen_rna

        gen_protein = options.get('gen_protein', True)
        assert(isinstance(gen_protein,bool))
        options['gen_protein'] = gen_protein

        gen_complexes = options.get('gen_complexes', True)
        assert(isinstance(gen_complexes,bool))
        options['gen_complexes'] = gen_complexes

        gen_concentrations = options.get('gen_concentrations', True)
        assert(isinstance(gen_concentrations,bool))
        options['gen_concentrations'] = gen_concentrations

        gen_observables = options.get('gen_observables', True)
        assert(isinstance(gen_observables,bool))
        options['gen_observables'] = gen_observables

        gen_kb_reactions = options.get('gen_kb_reactions', True)
        assert(isinstance(gen_kb_reactions,bool))
        options['gen_kb_reactions'] = gen_kb_reactions

        gen_kb_rate_laws = options.get('gen_kb_rate_laws', True)
        assert(isinstance(gen_kb_rate_laws,bool))
        options['gen_kb_rate_laws'] = gen_kb_rate_laws

    def gen_compartments(self):
        cell = self.knowledge_base.cell
        model = self.model

        # Create default compartments
        self.model.compartments.create(id='c', name='Cytosol', initial_volume=1E-15)
        self.model.compartments.create(id='e', name='Extracellular space', initial_volume=1E-10)

        # Generate the compartments that are defined in the knowledge base
        # TODO: currently no volume info in KB, talk to YH
        # for compartment_kb in self.knowledge_base.cell.compartments:
        #    compartment_model = self.model.compartments.get_or_create(id=kb_comp.id)
        #    compartment_model.name = compartment_kb.name
        #    compartment_model.initial_volume = compartment_kb.properties.get_one(id='volume').value

    def gen_parameters(self):
        # Create parameters
        self.model.parameters.create(id='fraction_dry_weight',
            value = self.knowledge_base.cell.properties.get_one(id='fraction_dry_weight').value)
        self.model.parameters.get_or_create(id='fraction_dry_weight').units = 'dimensionless'

        self.model.parameters.create(id='cell_cycle_length',
            value = self.knowledge_base.cell.properties.get_one(id='cell_cycle_length').value)
        self.model.parameters.get_or_create(id='cell_cycle_length').units = 's'

    def gen_metabolic_species(self):
        """ Generate all metabolic species in the cytosol """
        cytosol = self.model.compartments.get_one(id='c')
        metabolites = self.knowledge_base.cell.species_types.get(
            __type=wc_kb.core.MetaboliteSpeciesType)

        # get or create metabolite species
        for kb_met in metabolites:
            self.gen_a_specie(kb_met, cytosol)

    def gen_rna(self):
        '''Generate RNAs in wc_lang model from knowledge base '''
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]


        # get or create RNA species
        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna in rnas:
            species_type = model.species_types.get_or_create(id=rna.id)
            if not species_type.name:
                species_type.name = rna.name
                species_type.type = wc_lang.SpeciesTypeType.rna
                species_type.structure = rna.get_seq()
                species_type.empirical_formula = rna.get_empirical_formula()
                species_type.molecular_weight = rna.get_mol_wt()
                species_type.charge = rna.get_charge()
                species_type.comments = rna.comments
                species = species_type.species.get_or_create(compartment=cytosol)
                
    def gen_protein(self):
        '''Generate proteins in wc_lang model from knowledge base '''

        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType):

            species_type = self.model.species_types.get_or_create(
                id=protein.id)
            if not species_type.name:
                # Add functional form of protein
                species_type.name = protein.name
                species_type.type = wc_lang.SpeciesTypeType.protein
                species_type.structure = protein.get_seq()
                species_type.empirical_formula = protein.get_empirical_formula()
                species_type.molecular_weight = protein.get_mol_wt()
                species_type.charge = protein.get_charge()
                species_type.comments = protein.comments

                species = species_type.species.get_or_create(
                    compartment=cytosol)
                
    def gen_complexes(self):
        '''Generate complexes in wc_lang model from knowledge base '''
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

        for comp in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ComplexSpeciesType):
            species_type = model.species_types.get_or_create(id=comp.id)

            if not species_type.name:
                species_type.name = comp.name
                species_type.type = wc_lang.SpeciesTypeType.pseudo_species
                species_type.empirical_formula = comp.get_empirical_formula()
                species_type.molecular_weight = comp.get_mol_wt()
                species_type.charge = comp.get_charge()

                species = species_type.species.get_or_create(compartment=cytosol)
                
    def gen_concentrations(self):
        '''Generate concentrations in wc_lang model from knowledge base '''
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

        for conc in self.knowledge_base.cell.concentrations:
            species_type = model.species_types.get_or_create(id=conc.species.species_type.id)

            if not species_type.name:
                species_type.name = conc.species.species_type.name                
                species_type.comments = conc.species.species_type.comments
                species_type.references = conc.species.species_type.references
                
                species = species_type.species.get_or_create(compartment=cytosol)
                species.concentration = wc_lang.Concentration(
                    value=conc.value, units=wc_lang.ConcentrationUnit.M,
                    comments=conc.comments, references=conc.references) 

    def gen_observables(self):
        '''Generate observables in wc_lang model from knowledge base '''
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]
        observable_references = {Species:{}, Observable:{}}
        for kb_observable in self.knowledge_base.cell.observables:
            model_observable = self.model.observables.get_or_create(
                id=kb_observable.id)

            obs_expr_parts = []
            if not model_observable.name:
                model_observable.name = kb_observable.name
                for kb_species_coefficient in kb_observable.species:
                    kb_species = kb_species_coefficient.species
                    kb_species_type = kb_species.species_type
                    kb_compartment = kb_species.compartment
                    model_species_type = model.species_types.get_one(
                        id=kb_species_type.id)
                    model_species = model_species_type.species.get_one(
                        compartment=model.compartments.get_one(id=kb_compartment.id))
                    observable_references[Species][model_species.get_id()] = model_species
                    model_coefficient = kb_species_coefficient.coefficient
                    obs_expr_parts.append("{}*{}".format(model_coefficient, model_species.get_id()))

                for kb_observable_observable in kb_observable.observables:
                    model_observable_observable = model.observables.get_or_create(
                        id=kb_observable_observable.id)
                    obs_expr_parts.append("{}*{}".format(kb_observable_observable.coefficient, kb_observable_observable.id))
                    observable_references[Observable][model_observable_observable.id] = model_observable_observable
                obs_expr, e = ExpressionMethods.make_expression_obj(Observable,
                    ' + '.join(obs_expr_parts), observable_references)
                assert e is None, "cannot deserialize ObservableExpression: {}".format(e)
                model_observable.expression = obs_expr

    def gen_kb_reactions(self):
        """ Generate reactions encoded within KB """

        for kb_rxn in self.knowledge_base.cell.reactions:
            # if species are metabolites, create lang reaction in metabolism
            # todo: generalize to all submodels
            if self.get_species_type_types(kb_rxn) == set([wc_kb.core.MetaboliteSpeciesType]):
                lang_rxn = self.submodel.reactions.create(
                    id=kb_rxn.id,
                    name=kb_rxn.name,
                    reversible=kb_rxn.reversible,
                    comments=kb_rxn.comments)
                for participant in kb_rxn.participants:
                    kb_species = participant.species
                    lang_species_type = self.model.species_types.get_one(
                        id=kb_species.species_type.id)
                    lang_compartment = self.model.compartments.get_one(
                        id=kb_species.compartment.id)
                    lang_species = lang_species_type.species.get_one(
                        compartment=lang_compartment)
                    # ensure that species are present in extracellular space
                    if lang_species is None:
                        lang_species = self.gen_a_specie(kb_species.species_type, lang_compartment)
                    lang_rxn.participants.add(
                        lang_species.species_coefficients.get_or_create(
                            coefficient=participant.coefficient))

    def gen_kb_rate_laws(self):

        """ Generate rate laws for reactions encoded in KB """
        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='metabolism')
        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')

        for kb_rxn in self.knowledge_base.cell.reactions:
            if self.get_species_type_types(kb_rxn) == set([wc_kb.core.MetaboliteSpeciesType]):
                lang_rxn = self.submodel.reactions.get_one(id=kb_rxn.id)
                for kb_rate_law in kb_rxn.rate_laws:
                    lang_rate_law = wc_lang.RateLaw(k_cat=kb_rate_law.k_cat,
                        k_m=kb_rate_law.k_m,
                        comments=kb_rate_law.comments)
                    lang_rxn.rate_laws.add(lang_rate_law)
                    if kb_rate_law.direction == wc_kb.RateLawDirection.forward:
                        lang_rate_law.direction = wc_lang.RateLawDirection.forward
                    elif kb_rate_law.direction == wc_kb.RateLawDirection.backward:  # pragma branch not covered
                        lang_rate_law.direction = wc_lang.RateLawDirection.backward
                    lang_rate_law.equation = wc_lang.RateLawEquation(
                        expression=kb_rate_law.equation.expression
                    )

                    for kb_modifier in kb_rate_law.equation.modifiers:
                        lang_species_type = self.model.species_types.get_one(
                            id=kb_modifier.species_type.id)
                        lang_species = lang_species_type.species.get_one(
                            compartment=c)
                        lang_rate_law.equation.modifiers.add(lang_species)

    def gen_a_specie(self, kb_metabolite, lang_compartment):
        """ Generate a species in a particular compartment

        Args:
            kb_metabolite (:obj:`wc_kb.MetaboliteSpeciesType`): a knowledgebase metabolite
            lang_compartment (:obj:`wc_lang.Compartment`): the wc_lang compartment containing the species

        Returns:
            :obj:`wc_lang.Species`: the species that was found or created
        """
        species_type = self.model.species_types.get_or_create(id=kb_metabolite.id)
        species_type.name = kb_metabolite.name
        species_type.type = wc_lang.SpeciesTypeType.metabolite
        species_type.structure = kb_metabolite.structure
        species_type.empirical_formula = kb_metabolite.get_empirical_formula()
        species_type.molecular_weight = kb_metabolite.get_mol_wt()
        species_type.charge = kb_metabolite.get_charge()
        species_type.comments = kb_metabolite.comments
        species = species_type.species.get_or_create(
            compartment=lang_compartment)
        
        return species

    def get_species_type_types(self, kb_rxn):
        """ Obtain the species type types used by a kb reaction """
        species_type_types = set()
        for participant in kb_rxn.participants:
            kb_species_type = participant.species.species_type
            species_type_types.add(type(kb_species_type))
        return species_type_types
