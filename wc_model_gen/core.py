""" Base classes for generating :obj:`wc_lang`-formatted models from a knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import abc
import six
import wc_kb
import wc_lang
import math
import numpy
import wc_utils.util.string


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

        if rate_law_dynamics=='phenomenological':
            self.gen_phenom_rates()
        elif rate_law_dynamics=='mechanistic':
            self.gen_mechanistic_rates()
        else:
            raise Exception('Invalid rate law option selected.')

    def gen_phenom_rate_law_eq(self, specie_type_kb, reaction, half_life, cell_cycle_length):
        cytosol = self.model.compartments.get_one(id='c')
        specie_type_model = self.model.species_types.get_one(id=specie_type_kb.id)
        specie_model = specie_type_model.species.get_one(compartment=cytosol)

        rate_law = reaction.rate_laws.create()
        rate_law.direction = wc_lang.RateLawDirection.forward
        expression = '({} / {} + {} / {}) * {}'.format(numpy.log(2), half_life,
                                                       numpy.log(2), cell_cycle_length,
                                                       specie_model.id())

        rate_law.equation = wc_lang.RateLawEquation(expression = expression)
        rate_law.equation.modifiers.append(specie_model)

    def gen_mechanistic_rate_law_eq(self, submodel, specie_type_kb, reaction, beta, half_life, cell_cycle_length):

        cytosol_kb    = self.knowledge_base.cell.compartments.get_one(id='c')
        cytosol_model = self.model.compartments.get_one(id='c')
        specie_type_model = self.model.species_types.get_one(id=specie_type_kb.id)
        specie_model = specie_type_model.species.get_one(compartment=cytosol_model)

        expression = 'k_cat*'
        modifiers = []
        rate_avg = ''

        for participant in reaction.participants:
            if participant.coefficient < 0:
                avg_conc = (3/2)*participant.species.concentration.value
                modifiers.append(participant.species)
                rate_avg   += '({}/({}+({}*{})))*'.format(avg_conc, avg_conc, beta, avg_conc)
                expression += '({}/({}+({}*{})))*'.format(participant.species.id(),
                                                          participant.species.id(),
                                                          beta,
                                                          participant.species.concentration.value)

        # Clip off trailing '*' character
        expression = expression[:-1]
        rate_avg = rate_avg[:-1]

        # Check if RL eq already exists
        for sm_reaction in submodel.reactions:
            for rl in sm_reaction.rate_laws:
                if rl.equation.expression == expression:
                    rate_law_equation = rl.equation
                    break

        if 'rate_law_equation' not in locals():
            rate_law_equation = wc_lang.RateLawEquation(expression=expression, modifiers=modifiers)

        # Create rate law
        rate_law = reaction.rate_laws.create()
        rate_law.direction = wc_lang.RateLawDirection.forward
        rate_law.equation = rate_law_equation

        # Calculate k_cat
        exp_expression = '({}*(1/{}+1/{})*{})'.format(
                            numpy.log(2),
                            self.knowledge_base.cell.properties.get_one(id='cell_cycle_length').value,
                            half_life,
                            3/2*specie_type_kb.species.get_one(compartment=cytosol_kb).concentrations.value)

        rate_law.k_cat = eval(exp_expression) / eval(rate_avg)

    def calc_mean_half_life(self, species_types_kb):

        half_lifes=[]
        for species_type_kb in species_types_kb:
            if (isinstance(species_type_kb.half_life, float) and not species_type_kb.half_life==0 and not math.isnan(species_type_kb.half_life)):
                half_lifes.append(species_type_kb.half_life)
        avg_half_life = numpy.mean(half_lifes)

        return avg_half_life
