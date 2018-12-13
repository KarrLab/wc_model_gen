from wc_model_gen import prokaryote
import os
import shutil
import tempfile
import unittest
import wc_kb
import wc_kb_gen
import wc_lang


class ModelIOTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir)

    def test_generate_read_write(self):
        kb = wc_kb.io.Reader().run('tests/fixtures/test_broken_kb.xlsx',
                                   'tests/fixtures/test_broken_seq.fna',
                                   strict=False)

        """ Generate and write models from KBs """
        model = prokaryote.ProkaryoteModelGenerator(knowledge_base=kb).run()
        wc_lang.io.Writer().run(model, os.path.join(self.dir, 'model.xlsx'), set_repo_metadata_from_path=False)

        """ Read back KBs from disk """
        model_read = wc_lang.io.Reader().run(os.path.join(self.dir, 'model.xlsx'))
        self.assertTrue(model.is_equal(model_read))
