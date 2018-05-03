""" Tests of the API

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-04-19
:Copyright: 2018, Karr Lab
:License: MIT
"""

import types
import unittest
import wc_model_gen


class ApiTestCase(unittest.TestCase):
    def test(self):
        self.assertIsInstance(wc_model_gen.InitalizeModel, type)
