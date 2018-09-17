import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_lang
import wc_kb_gen
import wc_kb
import unittest
import tempfile
import shutil
import os

class KbAndModelIOTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_generate_read_write(self):

        """ Generate KBs """
        kb = wc_kb_gen.random.RandomKbGenerator(options={
                     'component': {
                         'GenomeGenerator': {
                             'genetic_code': 'reduced',
                             'mean_num_genes': 40,
                             'mean_gene_len': 70,
                             'num_tRNA': 4,
                             'num_rRNA': 0,
                             'num_ncRNA': 0,
                             'min_prots': 5,
                             'translation_table': 4,
                             'mean_copy_number': 100,
                             'mean_half_life': 100},
                         'PropertiesGenerator': {
                             'mean_cell_cycle_length': 100},
                         'ObservablesGenerator': {
                             'genetic_code': 'reduced'},
                         }}).run()

        self.assertIsInstance(kb, wc_kb.core.KnowledgeBase)

        """ Generate  and write models from KBs """
        model         = prokaryote.ProkaryoteModelGenerator(knowledge_base=kb).run()
        wc_lang.io.Writer().run(model, os.path.join(self.dir, 'model.xlsx'), set_repo_metadata_from_path=False)

        self.assertTrue(os.path.isfile(os.path.join(self.dir, 'model.xlsx')))

        """ Read back KBs from disk """
        #model_read = wc_lang.io.Reader().run('/media/sf_VM_share/model.xlsx')
