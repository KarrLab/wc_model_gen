import wc_model_gen.prokaryote as prokaryote
import unittest
import wc_lang
import wc_kb_gen
import wc_kb
import unittest
import tempfile
import shutil
import os

class ModelIOTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_generate_read_write(self):

        kb = wc_kb.io.Reader().run('tests/fixtures/test_broken.xlsx',
                                   'tests/fixtures/test_broken_seq.fna',
                                        strict=False)

        self.assertIsInstance(kb, wc_kb.core.KnowledgeBase)

        """ Generate  and write models from KBs """
        model         = prokaryote.ProkaryoteModelGenerator(knowledge_base=kb).run()
        wc_lang.io.Writer().run(model, os.path.join(self.dir, 'model.xlsx'), set_repo_metadata_from_path=False)
        self.assertIsInstance(model, wc_lang.core.Model)
        self.assertTrue(os.path.isfile(os.path.join(self.dir, 'model.xlsx')))

        """ Read back KBs from disk """
        model_read = wc_lang.io.Reader().run(os.path.join(self.dir, 'model.xlsx'))
        self.assertIsInstance(model_read, wc_lang.core.Model)
