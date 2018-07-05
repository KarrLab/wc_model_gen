""" Generating wc_lang formatted models from knowledge base.
:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
Model specification:
    Submodel is based on 'Physiology of the bacterial cell' by Neidhardt, Ingraham and Schaechter 1990.
    Translation is represented as 3 lumped reactions:
        Initation:   [c]: complex_70S_IA + GTP + IF1 + IF2 + IF3  ==> complex_70S_A + IF1 + IF2 + IF3 + GDP + P
        Elongation:  [c]: complex_70S_A + (n_i*tRNA_i) + n*(EFtu + EFts + EFg) + 2n*GTP ==> PROT_att + 2n*(GDP + P)
        Termination: [c]: (RF1 OR RF2) + RF3 + GTP  ==> PROT + GDP + P
    Component IDs:
        IF1: prot_MPN187; IF2: prot_MPN155; IF3: prot_MPN115
        EFtu: prot_MPN665; EFts: prot_MPN631; EFg:  prot_MPN227
        30S (prot_MPN): 164, 190, 225 ,189, 178, 622, 660, 174, 230, 169, 208, 541, 296, 171, 446, 182, 228, 226, 179, 616
        50S (prot_MPN): 220, 538, 219, 617, 175, 183, 172, 192, 181, 658, 168, 117, 325, 170, 167, 176,
                           327, 624, 173, 165, 360, 540, 471, 069, 682, 116, 188, 166, 177, 180, 539, 231
    Notes:
        complex_70S_IA = 50S+30S (Inactive form of ribosome)
        KB also has EF-P (MPN029) that does not appear in Neidhart et al., how to incorporate?
        PROT_att: Prot attached to ribosome
        RF1: UAG; RF2: UGA; either works for UAA
        Only RF1 found within DB (MPN361)!
    TODO:
        - add physiological k_cat & k_m values
        - paramters are fixed within each geneartion methods, add optional dict argument input
        - Make RF factors sensitive to termination codon
        - systematize initaiton of species, e.g. currently ribosome complexes are manually added @ translation.gen_species()
        - add method to display reactions participants in string format (serialize dnw)
"""

from .parameters import ParametersGenerator
from .compartments import CompartmentsGenerator
from .metabolites import MetaboliteSpeciesGenerator
from .transcription import TranscriptionSubmodelGenerator
from .translation import TranslationSubmodelGenerator
from .degradation import DegradationSubmodelGenerator
from mycoplasma_pneumoniae import config
import wc_kb
import wc_model_gen


class ModelGenerator(wc_model_gen.ModelGenerator):
    DEFAULT_COMPONENT_GENERATORS = (
        CompartmentsGenerator,
        MetaboliteSpeciesGenerator,
        ParametersGenerator,
    )

    def __init__(self, knowledge_base, component_generators=None, options=None):
        """
        Args:
            knowledge_base (:obj:`wc_kb.KnowledgeBase`): knowledge base
            component_generators (:obj:`list` of :obj:`wc_model_gen.ModelComponentGenerator`, optional): model component generators
            options (:obj:`dict`, optional): options
        """
        options = config.get_config({'model_gen': options or {}})['model_gen']
        super(ModelGenerator, self).__init__(knowledge_base, component_generators=component_generators, options=options)
