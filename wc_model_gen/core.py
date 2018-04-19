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
import wc_lang.core

class ModelGenerator(object):
    """ Generating a model instance (:obj:`wc_lang.core.Model`)

    Attributes:
        component_generators (:obj:`list` of :obj:`ModelComponentGenerator`): model component generators
        options (:obj:`dict`, optional): dictionary of options whose keys are methods and values are
            optional arguments to the methods
    """

    DEFAULT_MODEL_COMPONENTS = () # todo: run default components
    DEFAULT_MODEL_GENERATOR_VERSION = '0.0.1'
    DEFAULT_MODEL_VERSION = '0.0.1'
    DEFAULT_MODEL_ID = 'test_model'

    def __init__(self, knowledge_base, components=None, options=None, version=None):
        """
        Args:
            component_generators (:obj:`tuple` of :obj:`ModelComponentGenerator`): model component generators
            options (:obj:`dict`, optional): dictionary of options whose keys are method names and values are
                optional arguments to the methods
        """

        self.knowledge_base = knowledge_base
        self.components = components or self.DEFAULT_MODEL_COMPONENTS
        self.version = version or self.DEFAULT_MODEL_GENERATOR_VERSION
        self.options = options or {}

    def run(self, id=None, version=None):
        """ Generate a wc_lang model from a wc_kb knowledge base
        Args:
            id (:obj:`str`): model id

        Returns:
            :obj:`wc_lang.core.Model`: model
        """

        model = wc_lang.core.Model()
        model.id = id or self.DEFAULT_MODEL_ID
        model.version = version or self.DEFAULT_MODEL_VERSION

        return model


class ModelComponentGenerator(six.with_metaclass(abc.ABCMeta, object)):
    """ Base class for model component generator classes

    Attributes:
        knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
        model (:obj:`wc_lang.core.Model`): model
    """

    def __init__(self, knowledge_base, model):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            model (:obj:`wc_lang.core.Model`): model
        """
        self.knowledge_base = knowledge_base
        self.model = model

    @abc.abstractmethod
    def run(self):
        """ Generate species associated with submodel """
        pass  # pragma: no cover
