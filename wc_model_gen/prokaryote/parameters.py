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
        self.model.parameters.create(id='fractionDryWeight', value=0.7)
        self.model.parameters.create(id='cellCycleLength', value=100)
