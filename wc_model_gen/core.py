""" Base classes for generating wc_lang-formatted models from a knowledge base.

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


class ModelGenerator(object):
    """ Generating a model instance (:obj:`wc_lang.Model`)

    Attributes:
        component_generators (:obj:`list` of :obj:`ModelComponentGenerator`): model component generators
        options (:obj:`dict`, optional): dictionary of options whose keys are methods and values are
            optional arguments to the methods
    """

    DEFAULT_COMPONENTS = ()  # todo: run default components

    def __init__(self, knowledge_base, components=None, options=None):
        """
        Args:
            component_generators (:obj:`tuple` of :obj:`ModelComponentGenerator`): model component generators
            options (:obj:`dict`, optional): dictionary of options whose keys are method names and values are
                optional arguments to the methods
        """

        self.knowledge_base = knowledge_base
        self.components = components or self.DEFAULT_COMPONENTS
        self.options = options or {}

    def run(self, id=None, version=None):
        """ Generate a wc_lang model from a :obj:`wc_kb` knowledge base
        Args:
            id (:obj:`str`): model id

        Returns:
            :obj:`wc_lang.Model`: model
        """
        model = wc_lang.Model()
        model.id = self.options.get('id', None)
        model.version = self.options.get('version', None)

        # run component generators
        component_options = self.options.get('component', {})
        for component in self.components:
            options = component_options.get(component.__name__, {})
            component(self.knowledge_base, model, options=options).run()

        # return model
        return model


class ModelComponentGenerator(six.with_metaclass(abc.ABCMeta, object)):
    """ Base class for model component generator classes

    Attributes:
        knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
        model (:obj:`wc_lang.Model`): model
    """

    def __init__(self, knowledge_base, model, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            model (:obj:`wc_lang.Model`): model
            model (:obj:`dict`, optional): options
        """
        self.knowledge_base = knowledge_base
        self.model = model
        self.options = options

    @abc.abstractmethod
    def run(self):
        """ Generate species associated with submodel """
        pass  # pragma: no cover
