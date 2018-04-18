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
    """ Machinery for generating a model (an instance of :obj:`wc_lang.core.Model`) from a knowledge base
    (an instance of :obj:`wc_kb.KnowledgeBase`)

    Attributes:
        knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
        component_generators (:obj:`list` of :obj:`ModelComponentGenerator`): model component generators
        options (:obj:`dict`, optional): dictionary of options whose keys are methods and values are
            optional arguments to the methods
    """

    # Compartments + parameters are neccessary for a valid wc_lang model
    DEFAULT_MODEL_GENERATOR_VERSION = '0.0.1'
    DEFAULT_MODEL_COMPONENTS = () #(CompartmentsGenerator, ParametersGenerator)
    DEFAULT_MODEL_ID = 'test_model'

    def __init__(self, knowledge_base, components=None, options=None, version=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            component_generators (:obj:`tuple` of :obj:`ModelComponentGenerator`): model component generators
            options (:obj:`dict`, optional): dictionary of options whose keys are method names and values are
                optional arguments to the methods
        """

        self.knowledge_base = knowledge_base
        self.components = components or self.DEFAULT_MODEL_COMPONENTS
        self.version = version or self.DEFAULT_MODEL_GENERATOR_VERSION
        self.options = options or {}

    def run(self, id=None):
        """ Generate a wc_lang model from a wc_kb knowledge base
        Args:
            id (:obj:`str`): model id

        Returns:
            :obj:`wc_lang.core.Model`: model
        """

        model = wc_lang.core.Model()
        model.id = id or self.DEFAULT_MODEL_ID
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
        """ Run the generator """
        pass  # pragma: no cover
