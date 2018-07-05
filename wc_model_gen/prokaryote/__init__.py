import pkg_resources

with open(pkg_resources.resource_filename('wc_model_gen', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()

from .parameters import ParametersGenerator
from .compartments import CompartmentsGenerator
from .metabolites import MetaboliteSpeciesGenerator
from .transcription import TranscriptionSubmodelGenerator
from .translation import TranslationSubmodelGenerator
from .degradation import DegradationSubmodelGenerator
