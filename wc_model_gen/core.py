

""" Base classes for generating :obj:`wc_lang`-formatted models from a knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import abc
import math
import numpy
import os
import scipy.constants
import six
import time
import wc_kb
import wc_lang
import wc_utils.util.string
from wc_lang.util import get_model_summary
from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults
from wc_utils.util.ontology import wcm_ontology
from wc_utils.util.units import unit_registry


class ModelGenerator(object):
    """ Generator for models (:obj:`wc_lang.Model`)

    Attributes:
        knowledge_base (:obj:`wc_kb.core.KnowledgeBase`): knowledge base
        component_generators (:obj:`list` of :obj:`ModelComponentGenerator`): model component generators
        options (:obj:`dict`, optional): dictionary of options whose keys are the names of component
            generator classes and whose values are dictionaries of options for the component generator
            classes
    """

    DEFAULT_COMPONENT_GENERATORS = ()

    def __init__(self, knowledge_base, component_generators=None, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.core.KnowledgeBase`): knowledge base
            component_generators (:obj:`list` of :obj:`ModelComponentGenerator`, optional): model component generators
            options (:obj:`dict`, optional): dictionary of options whose keys are the names of component
                generator classes and whose values are dictionaries of options for the component generator
                classes
        """

        self.knowledge_base = knowledge_base

        if not component_generators:
            component_generators = list(self.DEFAULT_COMPONENT_GENERATORS)
        self.component_generators = component_generators

        self.options = options or {}
        self.clean_and_validate_options()

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        id = options.get('id', None)
        assert(isinstance(id, str) or id is None)
        options['id'] = id

        name = options.get('name', None)
        assert(isinstance(name, str) or name is None)
        options['name'] = name

        version = options.get('version', None)
        assert(isinstance(version, str) or version is None)
        options['version'] = version

    def run(self):
        """ Generate a :obj:`wc_lang` model from a :obj:`wc_kb` knowledge base

        Returns:
            :obj:`wc_lang.Model`: model
        """
        model = wc_lang.Model()
        model.id = self.options.get('id')
        model.name = self.options.get('name')
        model.version = self.options.get('version')

        component_options = self.options.get('component', {})
        for component_generator in self.component_generators:
            options = component_options.get(component_generator.__name__, {})
            component_generator(self.knowledge_base, model, options=options).run()

        return model

    """ Not sure what is the best place for the following static methods """
    @staticmethod
    def gen_rand_min_model_kb(name=None):
        """ Generates a random min model KB """

        kb = wc_kb_gen.random.RandomKbGenerator(options={
            'component': {
                'GenomeGenerator': {
                    'genetic_code': 'reduced',
                    'num_genes': 16,
                    'mean_gene_len': 50,
                    'num_tRNA': 4,
                    'num_rRNA': 0,
                    'num_ncRNA': 0,
                    'min_prots': 5,
                    'translation_table': 4,
                    'mean_rna_copy_number': 100,
                    'mean_protein_copy_number': 100},
                'PropertiesGenerator': {'mean_doubling_time': 100},
                'ObservablesGenerator': {'genetic_code': 'reduced'}}}).run()

        cell = kb.cell
        cytosol = cell.compartments.get_one(id='c')
        cell.species_types.get_one(id='amp').species.get_one(compartment=cytosol).concentration = \
            wc_kb.core.Concentration(cell=cell, value='0.000008333')
        cell.species_types.get_one(id='cmp').species.get_one(compartment=cytosol).concentration = \
            wc_kb.core.Concentration(cell=cell, value='0.000008333')
        cell.species_types.get_one(id='gmp').species.get_one(compartment=cytosol).concentration = \
            wc_kb.core.Concentration(cell=cell, value='0.000008333')
        cell.species_types.get_one(id='ump').species.get_one(compartment=cytosol).concentration = \
            wc_kb.core.Concentration(cell=cell, value='0.000008333')
        cell.species_types.get_one(id='atp').species.get_one(compartment=cytosol).concentration = \
            wc_kb.core.Concentration(cell=cell, value='0.000008333')
        cell.species_types.get_one(id='ctp').species.get_one(compartment=cytosol).concentration = \
            wc_kb.core.Concentration(cell=cell, value='0.000008333')
        cell.species_types.get_one(id='gtp').species.get_one(compartment=cytosol).concentration = \
            wc_kb.core.Concentration(cell=cell, value='0.000008333')
        cell.species_types.get_one(id='utp').species.get_one(compartment=cytosol).concentration = \
            wc_kb.core.Concentration(cell=cell, value='0.000008333')

        if name:
            core_path = '/media/sf_VM_share/kbs/' + str(name) + '.xlsx'
            seq_path = '/media/sf_VM_share/kbs/' + str(name) + '_seq.fna'

            wc_kb.io.Writer().run(core_path, kb,
                                  seq_path=seq_path,
                                  set_repo_metadata_from_path=False)

        return kb

    @staticmethod
    def run_model(model, results_dir, checkpoint_period=5, end_time=100):
        """ Simulates model """

        if not os.path.exists(results_dir):
            os.makedirs(results_dir)

        simulation = Simulation(model)
        results = simulation.run(end_time, results_dir, checkpoint_period)

        return results

    @staticmethod
    def analyze_model(self, results):
        """ Prints the standard analysis of simulation results """

        num_events = results[0]
        run_results_dir = results[1]
        run_results = RunResults(run_results_dir)

        rna_ids = []
        df = run_results.get('populations')
        for rna in model.species_types.get(type=wcm_ontology['WCM:RNA']): #RNA
            rna_ids.append(rna.species[0].id)

        txt = 'Init copy number mean={}; std={} \n'.format(round(np.mean(df.loc[0.0, rna_ids].values), 2),
                                                           round(np.std(df.loc[0.0, rna_ids].values), 2))
        txt += 'Final copy number mean={}; std={}'.format(round(np.mean(df.loc[100.0, rna_ids].values), 2),
                                                          round(np.std(df.loc[100.0, rna_ids].values), 2))
        print(txt)
        print(df['h[c]'], '\n')
        print(df[['amp[c]', 'cmp[c]', 'gmp[c]', 'ump[c]']], '\n')
        print(df[['atp[c]', 'ctp[c]', 'gtp[c]', 'utp[c]']], '\n')


class ModelComponentGenerator(six.with_metaclass(abc.ABCMeta, object)):
    """ Abstract base class for model component generators

    Attributes:
        knowledge_base (:obj:`wc_kb.core.KnowledgeBase`): knowledge base
        model (:obj:`wc_lang.Model`): model
        options (:obj:`dict`, optional): options
    """

    def __init__(self, knowledge_base, model, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.core.KnowledgeBase`): knowledge base
            model (:obj:`wc_lang.Model`): model
            options (:obj:`dict`, optional): options
        """
        self.knowledge_base = knowledge_base
        self.model = model
        self.options = options or {}

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        pass  # pragma: no cover

    @abc.abstractmethod
    def run(self):
        """ Generate model components """
        pass  # pragma: no cover


class SubmodelGenerator(ModelComponentGenerator):
    """ Base class for submodel generators

    Attributes:
        knowledge_base (:obj:`wc_kb.core.KnowledgeBase`): knowledge base
        model (:obj:`wc_lang.Model`): model
        submodel (:obj:`wc_lang.Submodel`): submodel
        options (:obj:`dict`, optional): options
    """

    def __init__(self, knowledge_base, model, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.core.KnowledgeBase`): knowledge base
            model (:obj:`wc_lang.Model`): model
            options (:obj:`dict`, optional): options
        """
        self.knowledge_base = knowledge_base
        self.model = model
        self.submodel = self.model.submodels.get_or_create(id=wc_utils.util.string.camel_case_to_snake_case(
            self.__class__.__name__.replace('SubmodelGenerator', '')))

        self.options = options or {}
        self.clean_and_validate_options()

        # Calculate numbers needed for model construction

    def run(self):
        """ Generate model components """
        self.clean_and_validate_options()
        self.gen_reactions()
        self.gen_rate_laws()

    def clean_and_validate_options(self):
        """ Apply default options and validate options """

        options = self.options

        rate_law_dynamics = options.get('rate_dynamics', 'phenomenological')
        assert(rate_law_dynamics in ['phenomenological', 'mechanistic'])
        options['rate_dynamics'] = rate_law_dynamics

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        pass  # pragma: no cover

    def gen_rate_laws(self):
        """ Choose which rate_law to generate """

        rate_law_dynamics = self.options.get('rate_dynamics')

        if rate_law_dynamics == 'phenomenological':
            self.gen_phenom_rates()
        elif rate_law_dynamics == 'mechanistic':
            self.gen_mechanistic_rates()
        else:
            raise Exception('Invalid rate law option.')

    def gen_phenom_rate_law_eq(self, specie_type_kb, reaction, half_life, mean_doubling_time):
        if reaction.id.endswith('_kb'):
            return

        model = self.model
        cytosol = model.compartments.get_one(id='c')
        species_type_model = model.species_types.get_one(id=specie_type_kb.id)
        species_model = species_type_model.species.get_one(compartment=cytosol)

        rate_law_model = model.rate_laws.create(
            reaction=reaction,
            direction=wc_lang.RateLawDirection.forward)
        rate_law_model.id = rate_law_model.gen_id()

        half_life_model = model.parameters.get_or_create(
            id='half_life_{}'.format(species_type_model.id),
            type=None,
            value=half_life,
            units=unit_registry.parse_units('s'))
        mean_doubling_time_model = model.parameters.get_or_create(
            id='mean_doubling_time',
            type=None,
            value=mean_doubling_time,
            units=unit_registry.parse_units('s'))
        molecule_units = model.parameters.get_or_create(
            id='molecule_units',
            type=None,
            value=1.,
            units=unit_registry.parse_units('molecule'))

        expression = '(log(2) / {} + log(2) / {}) / {} * {}'.format(half_life_model.id,
                                                                    mean_doubling_time_model.id,
                                                                    molecule_units.id,
                                                                    species_model.id)
        rate_law_model.expression, error = wc_lang.RateLawExpression.deserialize(expression, {
            wc_lang.Parameter: {
                half_life_model.id: half_life_model,
                mean_doubling_time_model.id: mean_doubling_time_model,
                molecule_units.id: molecule_units,
            },
            wc_lang.Species: {
                species_model.id: species_model,
            },
        })
        assert error is None, str(error)

    def gen_mechanistic_rate_law_eq(self, submodel, specie_type_kb, reaction, beta, half_life, mean_doubling_time):
        if reaction.id.endswith('_kb'):
            return

        kb = self.knowledge_base
        model = self.model
        cytosol_kb = kb.cell.compartments.get_one(id='c')
        cytosol_model = model.compartments.get_one(id='c')
        specie_type_model = model.species_types.get_one(id=specie_type_kb.id)
        specie_model = specie_type_model.species.get_one(compartment=cytosol_model)

        # Create rate law
        rate_law = model.rate_laws.create(
            reaction=reaction,
            direction=wc_lang.RateLawDirection.forward)
        rate_law.id = rate_law.gen_id()

        objects = {
            wc_lang.Species: {},
            wc_lang.Function: {},
            wc_lang.Parameter: {},
        }

        Avogadro = model.parameters.get_or_create(id='Avogadro',
                                                  type=None,
                                                     value=scipy.constants.Avogadro,
                                                     units=unit_registry.parse_units('molecule mol^-1'))
        objects[wc_lang.Parameter][Avogadro.id] = Avogadro

        model_k_cat = model.parameters.create(id='k_cat_{}'.format(reaction.id),
                                              type=wcm_ontology['WCM:k_cat'],
                                              value=1.,
                                              units=unit_registry.parse_units('s^-1'))
        objects[wc_lang.Parameter][model_k_cat.id] = model_k_cat

        expression_terms = []
        init_species_counts = {}

        init_compartment_masses = {}
        for comp in model.compartments:
            init_comp_mass = 0
            for species in comp.species:
                if species.distribution_init_concentration:
                    if species.distribution_init_concentration.units == unit_registry.parse_units('molecule'):
                        init_comp_mass += species.distribution_init_concentration.mean \
                            / scipy.constants.Avogadro \
                            * species.species_type.molecular_weight
                    elif species.distribution_init_concentration.units == unit_registry.parse_units('M'):
                        init_comp_mass += species.distribution_init_concentration.mean \
                            * comp.mean_init_volume \
                            * species.species_type.molecular_weight
                    else:
                        raise Exception('Unsupported concentration unit {}'.format(
                            species.distribution_init_concentration.units))
            init_compartment_masses[comp.id] = init_comp_mass

        for species in reaction.get_reactants():
            objects[wc_lang.Species][species.id] = species

            model_k_m = model.parameters.create(id='K_m_{}_{}'.format(reaction.id, species.species_type.id),
                                                type=wcm_ontology['WCM:K_m'],
                                                value=beta * species.distribution_init_concentration.mean,
                                                units=unit_registry.parse_units('M'))
            objects[wc_lang.Parameter][model_k_m.id] = model_k_m

            volume = species.compartment.init_density.function_expressions[0].function
            objects[wc_lang.Function][volume.id] = volume

            expression_terms.append('({} / ({} + {} * {} * {}))'.format(species.id,
                                                                        species.id,
                                                                        model_k_m.id, Avogadro.id,
                                                                        volume.id))

            # Calculate Avg copy number / concentration
            if species.distribution_init_concentration.units == unit_registry.parse_units('molecule'):
                init_species_counts[species.id] = species.distribution_init_concentration.mean
            elif species.distribution_init_concentration.units == unit_registry.parse_units('M'):
                init_species_counts[species.id] = species.distribution_init_concentration.mean \
                    * species.compartment.mean_init_volume \
                    * scipy.constants.Avogadro
            else:
                raise Exception('Unsupported units: {}'.format(species.distribution_init_concentration.units.name))

        expression = '{} * {}'.format(model_k_cat.id, ' * '.join(expression_terms))
        rate_law.expression, error = wc_lang.RateLawExpression.deserialize(expression, objects)
        assert error is None, str(error)

        # Calculate k_cat value
        init_species_conc = specie_type_kb.species.get_one(compartment=cytosol_kb).concentration.value
        int_species_cn = init_species_conc * cytosol_model.mean_init_volume * scipy.constants.Avogadro

        avg_deg_rate = math.log(2) * (1. / kb.cell.properties.get_one(id='mean_doubling_time').value + 1. / half_life) * int_species_cn
        rate_avg = rate_law.expression._parsed_expression.eval({
            wc_lang.Species: init_species_counts,
            wc_lang.Compartment: init_compartment_masses,
        })
        model_k_cat.value = avg_deg_rate / rate_avg
