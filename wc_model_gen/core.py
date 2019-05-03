""" Base classes for generating :obj:`wc_lang`-formatted models from a knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import abc
import os
import six
import time
import wc_kb
import wc_lang
import wc_utils.util.string
from wc_lang.util import get_model_summary
from wc_sim.multialgorithm.simulation import Simulation
from wc_sim.multialgorithm.run_results import RunResults
from wc_onto import onto as wc_ontology


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
            core_path = 'wc_model_gen/' + str(name) + '.xlsx'
            seq_path = 'wc_model_gen/' + str(name) + '_seq.fna'

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

        for rna in model.species_types.get(type=wc_ontology['WC:RNA']): #RNA
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
        self.calibrate_submodel()

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        pass  # pragma: no cover

    def gen_reactions(self):
        """ Generate reactions associated with the submodel """
        pass  # pragma: no cover

    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """
        pass  # pragma: no cover

    def calibrate_submodel(self):
        """ Calibrate the submodel using data in the KB """
        pass # pragma: no cover
