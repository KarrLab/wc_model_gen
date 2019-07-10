""" Generator for models based on KBs

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-07
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .initialize_model import InitializeModel
from .transcription import TranscriptionSubmodelGenerator
import wc_model_gen


class EukaryoteModelGenerator(wc_model_gen.ModelGenerator):
    """ Generator for submodels based on KBs

    Options:
    * id
    * name
    * version
    * component

        * InitializeModel
        * TranscriptionSubmodelGenerator
    """

    DEFAULT_COMPONENT_GENERATORS = (
        InitializeModel,
        TranscriptionSubmodelGenerator,
    )

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        id = options.get('id', 'eukaryote')
        assert(isinstance(id, str) or id is None)
        options['id'] = id

        name = options.get('name', 'Whole-cell model of eukaryote cells')
        assert(isinstance(name, str) or name is None)
        options['name'] = name

        version = options.get('version', wc_model_gen.__version__)
        assert(isinstance(version, str) or version is None)
        options['version'] = version
