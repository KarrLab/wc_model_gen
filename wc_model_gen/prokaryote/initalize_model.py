""" Initalize the construction of wc_lang-encoded models from wc_kb-encoded knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT

TODO:
- read cytosol volume from DB; currently there is only fractional volume?!
"""

from wc_lang import Species, Observable, ExpressionMethods
import numpy as np
import wc_model_gen
import wc_lang
import wc_kb
import math

class InitalizeModel(wc_model_gen.ModelComponentGenerator):
    """ Generate compartments """

    def run(self):
        self.gen_compartments()
        self.gen_parameters()
        self.clean_and_validate_options()
        options = self.options

        if options['gen_metabolites']:
            self.gen_metabolites()

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

        gen_metabolites = options.get('gen_metabolites', True)
        assert(isinstance(gen_metabolites, bool))
        options['gen_metabolites'] = gen_metabolites

        gen_rna = options.get('gen_rna', True)
        assert(isinstance(gen_rna, bool))
        options['gen_rna'] = gen_rna

        gen_protein = options.get('gen_protein', True)
        assert(isinstance(gen_protein, bool))
        options['gen_protein'] = gen_protein

        gen_complexes = options.get('gen_complexes', True)
        assert(isinstance(gen_complexes, bool))
        options['gen_complexes'] = gen_complexes

        gen_concentrations = options.get('gen_concentrations', True)
        assert(isinstance(gen_concentrations, bool))
        options['gen_concentrations'] = gen_concentrations

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

        # Create default compartments
        # TODO: currently no volume info in KB, talk to YH

        model.compartments.create(id='c', name='Cytosol', initial_volume=1E-15)
        #model.compartments.create(id='m', name='Cell membrane', initial_volume=1E-10)
        model.compartments.create(id='e', name='Extracellular space', initial_volume=1E-10)

    def gen_parameters(self):
        kb = self.knowledge_base
        model = self.model

        # Create parameters
        model.parameters.create(id='fraction_dry_weight',
                                value=kb.cell.properties.get_one(id='fraction_dry_weight').value)
        model.parameters.get_or_create(id='fraction_dry_weight').units = 'dimensionless'

        model.parameters.create(id='cell_cycle_length',
                                value=kb.cell.properties.get_one(id='cell_cycle_length').value)
        model.parameters.get_or_create(id='cell_cycle_length').units = 's'

        model.parameters.create(id='fractionDryWeight',
                                value=kb.cell.properties.get_one(id='fraction_dry_weight').value)
        model.parameters.get_or_create(id='fractionDryWeight').units = 'dimensionless'

        model.parameters.create(id='cellCycleLength',
                                value=kb.cell.properties.get_one(id='cell_cycle_length').value)
        model.parameters.get_or_create(id='cellCycleLength').units = 's'

    def gen_metabolites(self):
        """ Generate all metabolic species in the cytosol """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(
            __type=wc_kb.core.MetaboliteSpeciesType)

        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type)

    def gen_rna(self):
        """ Generate RNAs in model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.RnaSpeciesType)

        # Calculate average RNA half-life, will be used if value is missing
        half_lifes = []
        for kb_species_type in kb_species_types:
            if (isinstance(kb_species_type.half_life, float) and
                not kb_species_type.half_life == 0 and
                    not math.isnan(kb_species_type.half_life)):
                half_lifes.append(kb_species_type.half_life)

        avg_rna_half_life = np.mean(half_lifes)

        # Create RNA species
        for kb_species_type in kb_species_types:
            if math.isnan(kb_species_type.half_life) or kb_species_type.half_life == 0:
                # TODO: remove; KB should be treated as read only
                kb_species_type.half_life = avg_rna_half_life

            self.gen_species_type(kb_species_type, ['c'])

    def gen_protein(self):
        """ Generate proteins in model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)

        # Calculate average protein half-life, will be used if value is missing
        half_lifes = []
        for kb_species_type in kb_species_types:
            if isinstance(kb_species_type.half_life, float) and \
                    not kb_species_type.half_life == 0 and \
                    not math.isnan(kb_species_type.half_life):
                half_lifes.append(kb_species_type.half_life)

        avg_protein_half_life = np.mean(half_lifes)

        for kb_species_type in kb_species_types:
            if math.isnan(kb_species_type.half_life) or kb_species_type.half_life == 0:
                # TODO: remove; KB should be treated as read only
                kb_species_type.half_life = avg_protein_half_life

            self.gen_species_type(kb_species_type, ['c'])

    def gen_complexes(self):
        """ Generate complexes in model from knowledge base """
        kb = self.knowledge_base
        model = self.model

        kb_species_types = kb.cell.species_types.get(__type=wc_kb.core.ComplexSpeciesType)
        for kb_species_type in kb_species_types:
            self.gen_species_type(kb_species_type, ['c'])

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

        if isinstance(kb_species_type, wc_kb.core.MetaboliteSpeciesType):
            model_species_type.type = wc_lang.SpeciesTypeType.metabolite
            model_species_type.structure = kb_species_type.structure

        elif isinstance(kb_species_type, wc_kb.core.DnaSpeciesType):
            model_species_type.type = wc_lang.SpeciesTypeType.dna
            model_species_type.structure = kb_species_type.get_seq()

        elif isinstance(kb_species_type, wc_kb.prokaryote_schema.RnaSpeciesType):
            model_species_type.type = wc_lang.SpeciesTypeType.rna
            model_species_type.structure = kb_species_type.get_seq()

        elif isinstance(kb_species_type, wc_kb.prokaryote_schema.ProteinSpeciesType):
            model_species_type.type = wc_lang.SpeciesTypeType.protein
            model_species_type.structure = kb_species_type.get_seq()

        elif isinstance(kb_species_type, wc_kb.core.ComplexSpeciesType):
            model_species_type.type = wc_lang.SpeciesTypeType.protein
            model_species_type.structure = None

        else:
            raise ValueError('Unsupported species type: {}'.format(
                kb_species_type.__class__.__name__))

        model_species_type.empirical_formula = kb_species_type.get_empirical_formula()
        model_species_type.molecular_weight = kb_species_type.get_mol_wt()
        model_species_type.charge = kb_species_type.get_charge()
        model_species_type.comments = kb_species_type.comments
        compartment_ids = set([s.compartment.id for s in kb_species_type.species] +
                              (extra_compartment_ids or []))

        for compartment_id in compartment_ids:
            model_compartment = model.compartments.get_one(id=compartment_id)
            model_species = model_species_type.species.get_or_create(compartment=model_compartment)
            model_species.id = model_species.gen_id(model_species_type.id,
                                                    model_compartment.id)

        return model_species_type

    def gen_concentrations(self):
        """ Generate concentrations in model from knowledge base """
        kb = self.knowledge_base
        model = self.model
        cytosol = model.compartments.get_one(id='c')

        for conc in kb.cell.concentrations:
            species_comp_model = model.compartments.get_one(id=conc.species.compartment.id)

            species_type = model.species_types.get_or_create(id=conc.species.species_type.id)
            species = species_type.species.get_or_create(compartment=species_comp_model)
            species.id = species.gen_id(species.species_type.id, species.compartment.id)

            species.concentration = wc_lang.Concentration(
                value=conc.value, units=wc_lang.ConcentrationUnit.M,
                comments=conc.comments, references=conc.references)

    def gen_observables(self):
        """ Generate observables in model from knowledge base """
        kb = self.knowledge_base
        model = self.model
        cytosol = model.compartments.get(id='c')[0]
        observable_references = {Species: {}, Observable: {}}
        for kb_observable in kb.cell.observables:
            model_observable = model.observables.get_or_create(id=kb_observable.id)
            obs_expr_parts = []

            model_observable.name = kb_observable.name
            for kb_species_coefficient in kb_observable.species:
                kb_species = kb_species_coefficient.species
                kb_species_type = kb_species.species_type
                kb_compartment = kb_species.compartment
                model_species_type = model.species_types.get_one(
                    id=kb_species_type.id)
                model_species = model_species_type.species.get_one(
                    compartment=model.compartments.get_one(id=kb_compartment.id))
                observable_references[Species][model_species.id] = model_species
                model_coefficient = kb_species_coefficient.coefficient
                obs_expr_parts.append("{}*{}".format(model_coefficient, model_species.id))

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
        kb = self.knowledge_base
        model = self.model

        for kb_rxn in kb.cell.reactions:
            submodel_id = kb_rxn.submodel
            submodel = model.submodels.get_or_create(id=submodel_id)

            model_rxn = submodel.reactions.create(
                id=kb_rxn.id+'_fromKB',
                name=kb_rxn.name,
                reversible=kb_rxn.reversible,
                comments=kb_rxn.comments)

            for participant in kb_rxn.participants:
                kb_species = participant.species
                model_species_type = model.species_types.get_one(id=kb_species.species_type.id)
                model_compartment = model.compartments.get_one(id=kb_species.compartment.id)
                model_species = model_species_type.species.get_one(compartment=model_compartment)

                # ensure that species are present in extracellular space
                if model_species is None:
                    model_species = self.gen_species_type(kb_species.species_type, model_compartment)
                model_rxn.participants.add(
                    model_species.species_coefficients.get_or_create(coefficient=participant.coefficient))

    def gen_kb_rate_laws(self):
        """ Generate rate laws for reactions encoded in KB """
        kb = self.knowledge_base
        model = self.model

        c = model.compartments.get_one(id='c')
        e = model.compartments.get_one(id='e')

        for kb_rxn in kb.cell.reactions:
            submodel_id = kb_rxn.submodel
            submodel = model.submodels.get_or_create(id=submodel_id)
            model_rxn = submodel.reactions.get_one(id=kb_rxn.id+'_fromKB')

            for kb_rate_law in kb_rxn.rate_laws:

                model_rate_law = wc_lang.RateLaw(
                    k_cat=kb_rate_law.k_cat,
                    k_m=kb_rate_law.k_m,
                    comments=kb_rate_law.comments)

                model_rxn.rate_laws.add(model_rate_law)

                if kb_rate_law.direction == wc_kb.RateLawDirection.forward:
                    model_rate_law.direction = wc_lang.RateLawDirection.forward
                elif kb_rate_law.direction == wc_kb.RateLawDirection.backward:
                    model_rate_law.direction = wc_lang.RateLawDirection.backward

                model_rate_law.equation = wc_lang.RateLawEquation(expression=kb_rate_law.equation.expression)

                for kb_modifier in kb_rate_law.equation.modifiers:
                    model_species_type = model.species_types.get_one(id=kb_modifier.species_type.id)
                    model_species = model_species_type.species.get_one(compartment=c)
                    model_rate_law.equation.modifiers.add(model_species)
