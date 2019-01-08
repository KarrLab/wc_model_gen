from wc_model_gen import prokaryote
import os
import shutil
import tempfile
import unittest
import wc_kb
import wc_kb_gen
import wc_lang

kb = wc_kb.io.Reader().run('/home/balazs/Desktop/wc_model_generator/tests/fixtures/min_model_kb.xlsx',
                           '/home/balazs/Desktop/wc_model_generator/tests/fixtures/min_model_kb_seq.fna')

kb = kb[wc_kb.core.KnowledgeBase][0]

model = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=kb,
            component_generators=[prokaryote.InitalizeModel,
                                  prokaryote.TranscriptionSubmodelGenerator,
                                  prokaryote.RnaDegradationSubmodelGenerator,
                                  prokaryote.MetabolismSubmodelGenerator]).run()

wc_lang.io.Writer().run('/media/sf_VM_share/model_inspect.xlsx', model, set_repo_metadata_from_path=False)
