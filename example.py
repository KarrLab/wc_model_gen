from wc_model_gen import prokaryote
import wc_model_gen
import wc_lang
import wc_kb_gen
import wc_kb

rand_kb = wc_kb_gen.random.RandomKbGenerator(options={
             'component': {
                 'GenomeGenerator': {
                     'num_genes': 40,
                     'mean_gene_len': 50,
                     'num_ncRNA': 3,
                     'translation_table': 4,
                     'mean_copy_number': 2000,
                     'mean_half_life': 100
                 },
                 'PropertiesGenerator': {
                     'mean_cell_cycle_length': 300,
                 },
             },
         }).run()

model = prokaryote.ProkaryoteModelGenerator(
            knowledge_base=rand_kb,
            component_generators=[prokaryote.CompartmentsGenerator,
                                  prokaryote.ParametersGenerator,
                                  prokaryote.MetabolismSubmodelGenerator,
                                  prokaryote.TranscriptionSubmodelGenerator,
                                  prokaryote.RnaDegradationSubmodelGenerator]).run()

model.id = 'rand_kb'
model.version = '0.0.1'

wc_lang.io.Writer().run(model,'/media/sf_VM_share/models/rand_model_test.xlsx', set_repo_metadata_from_path=False)
