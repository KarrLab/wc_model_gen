""" Generating wc_lang formatted models from knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen


class ParametersGenerator(wc_model_gen.ModelComponentGenerator):
    """ Generate parameters """

    def run(self):
        self.model.parameters.create(id='fraction_dry_weight', value=0.7)
        self.model.parameters.create(id='cell_cycle_length', value=8 * 60 * 60)
