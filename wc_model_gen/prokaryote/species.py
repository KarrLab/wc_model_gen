""" Generate wc_lang species and wc_lang observables from the provided KB

:Author: Bilal Shaikh
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-08-01
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen
import wc_kb
import wc_lang
from wc_lang import Species, Observable, ExpressionMethods


class SpeciesGenerator(wc_model_gen.ModelComponentGenerator):

    def run(self):
        '''Generates model species from given KB'''
        self.gen_rna()
        self.gen_protein()
        self.gen_complexes()
        self.gen_observables()

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
                species = species_type.species.get_or_create(
                    compartment=cytosol)
                species.concentration = wc_lang.Concentration(
                    value=rna.concentration, units=wc_lang.ConcentrationUnit.M)

    def gen_protein(self):
        '''Generate proteins in wc_lang model from knowledge base '''

        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):

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

                species.concentration = wc_lang.Concentration(
                    value=protein.concentration, units=wc_lang.ConcentrationUnit.M)

    def gen_complexes(self):
        '''Generate complexes in wc_lang model from knowledge base '''
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]
        for comp in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ComplexSpeciesType):
            species_type = self.model.species_types.get_or_create(
                id=comp.id)
            if not species_type.name:
                species_type.name = comp.name
                species_type.type = wc_lang.SpeciesTypeType.pseudo_species
                species = species_type.species.get_or_create(
                    compartment=cytosol)

                species.concentration = wc_lang.Concentration(
                    value=comp.concentration, units=wc_lang.ConcentrationUnit.M)

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
