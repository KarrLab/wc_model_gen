""" Generating wc_lang formatted models from knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_model_gen


class CompartmentsGenerator(wc_model_gen.ModelComponentGenerator):
    """ Generate compartments """

    def run(self):
        self.model.compartments.create(id='c', name='Cytosol', initial_volume=4.60E-17)
        self.model.compartments.create(id='e', name='Extracellular space', initial_volume=1E-12)
