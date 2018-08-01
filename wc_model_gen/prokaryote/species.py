""" Generate wc_lang species and wc_lang observables from the provided KB

:Author: Bilal Shaikh
:Author Ashwin 
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen


class SpeciesGenerator(wc_model_gen.ModelComponentGenerator):
    def run(self):
        self.gen_rna()
        self.gen_protein()
        self.gen_observables()
        self.gen_complexes()

    def gen_rna(self):
        pass

    def gen_protein(self):
        pass

    def gen_complexes(self):
        pass

    def gen_observables(self):
        pass
