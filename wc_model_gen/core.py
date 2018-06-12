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
import wc_utils.util.string


class ModelGenerator(object):
    """ Generator for models (:obj:`wc_lang.Model`)

    Attributes:
        knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
        component_generators (:obj:`list` of :obj:`ModelComponentGenerator`): model component generators
        options (:obj:`dict`, optional): dictionary of options whose keys are the names of component
            generator classes and whose values are dictionaries of options for the component generator
            classes
    """

    DEFAULT_COMPONENT_GENERATORS = ()

    def __init__(self, knowledge_base, component_generators=None, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
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

    def run(self):
        """ Generate a :obj:`wc_lang` model from a :obj:`wc_kb` knowledge base

        Returns:
            :obj:`wc_lang.Model`: model
        """

        model = wc_lang.Model()
        model.id = self.options.get('id', None)
        model.version = self.options.get('version', None)

        component_options = self.options.get('component', {})
        for component_generator in self.component_generators:
            options = component_options.get(component_generator.__name__, {})
            component_generator(self.knowledge_base, model, options=options).run()

        return model


class ModelComponentGenerator(six.with_metaclass(abc.ABCMeta, object)):
    """ Base class for model component generators

    Attributes:
        knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
        model (:obj:`wc_lang.Model`): model
        options (:obj:`dict`, optional): options
    """

    def __init__(self, knowledge_base, model, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            model (:obj:`wc_lang.Model`): model
            options (:obj:`dict`, optional): options
        """
        self.knowledge_base = knowledge_base
        self.model = model
        self.options = options

    @abc.abstractmethod
    def run(self):
        """ Generate model components """
        pass  # pragma: no cover

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        pass  # pragma: no cover


class SubmodelGenerator(ModelComponentGenerator):
    """ Base class for submodel generators

    Attributes:
        knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
        model (:obj:`wc_lang.Model`): model
        submodel (:obj:`wc_lang.Submodel`): submodel
        options (:obj:`dict`, optional): options
    """

    def __init__(self, knowledge_base, model, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            model (:obj:`wc_lang.Model`): model            
            options (:obj:`dict`, optional): options
        """
        self.knowledge_base = knowledge_base
        self.model = model
        self.submodel = self.model.submodels.get_or_create(id=wc_utils.util.string.camel_case_to_snake_case(
            self.__class__.__name__.replace('SubmodelGenerator', '')))
        self.options = options

    def run(self):
        """ Generate model components """
        self.clean_and_validate_options()
        self.gen_compartments()
        self.gen_species()
        self.gen_reactions()
        self.gen_parameters()
        self.gen_rate_laws()

    def gen_compartments(self):
        """ Generate compartments associated with submodel """
        pass  # pragma: no cover

    def gen_species(self):
        """ Generate species associated with submodel """
        pass  # pragma: no cover

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        pass  # pragma: no cover

    def gen_parameters(self):
        """ Generate parameters associated with submodel """
        pass  # pragma: no cover

    def gen_rate_laws(self):
        """ Generate rate laws associated with submodel """
        pass  # pragma: no cover
