""" This file demonstrates how to use wc_model_gen to construct a wc_lang model from
an existing KB

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2019-03-15
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_model_gen import prokaryote
import wc_model_gen
import wc_lang
import wc_kb_gen
import wc_kb

# Read existing KB
kb = wc_kb.io.Reader().run('wc_model_gen/tests/fixtures/min_model_kb.xlsx',
                           'wc_model_gen/tests/fixtures/min_model_kb_seq.fna')
kb = kb[wc_kb.core.KnowledgeBase][0]


# Generate model, select which KB to use
model = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=kb,
            component_generators=[prokaryote.InitalizeModel,
                                  prokaryote.TranscriptionSubmodelGenerator,
                                  prokaryote.RnaDegradationSubmodelGenerator,
                                  prokaryote.MetabolismSubmodelGenerator]).run()

# Optionally save the resulting model as an excel file for inspection
#wc_lang.io.Writer().run('gen_model.xlsx', model, set_repo_metadata_from_path=False)
