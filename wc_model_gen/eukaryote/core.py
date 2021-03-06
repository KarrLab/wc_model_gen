""" Generator for models based on KBs

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2019-01-07
:Copyright: 2019, Karr Lab
:License: MIT
"""

from .complexation import ComplexationSubmodelGenerator
from .initialize_model import InitializeModel
from .metabolism import MetabolismSubmodelGenerator
from .protein_degradation import ProteinDegradationSubmodelGenerator
from .rna_degradation import RnaDegradationSubmodelGenerator
from .transcription import TranscriptionSubmodelGenerator
from .translation_translocation import TranslationTranslocationSubmodelGenerator
import wc_model_gen


class EukaryoteModelGenerator(wc_model_gen.ModelGenerator):
    """ Generator for submodels based on KBs

    Options:
    * id
    * name
    * version
    * component

        * InitializeModel
        * ComplexationSubmodelGenerator,
        * TranscriptionSubmodelGenerator,
        * RnaDegradationSubmodelGenerator
    """

    DEFAULT_COMPONENT_GENERATORS = (
        InitializeModel,
        ComplexationSubmodelGenerator,
        TranscriptionSubmodelGenerator,
        RnaDegradationSubmodelGenerator,
        TranslationTranslocationSubmodelGenerator,
        ProteinDegradationSubmodelGenerator,
        MetabolismSubmodelGenerator,
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
