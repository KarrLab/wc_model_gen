""" Generating wc_lang formatted models from knowledge base.

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT
"""

import wc_lang
import wc_model_gen


class MetaboliteSpeciesGenerator(wc_model_gen.ModelComponentGenerator):
    """ Generate species and cytosolic concentrations for each metabolite """

    def run(self):
        compartment = self.model.compartments.get_one(id='c')

        species_types = [
            ('ATP', 'C10H12N5O13P3', -4, 503.149),
            ('CTP', 'C9H12N3O14P3', -4, 479.123),
            ('GTP', 'C10H12N5O14P3', -4, 519.148),
            ('UTP', 'C9H11N2O15P3', -4, 480.107),
            ('ADP', 'C10H12N5O10P2', -3, 424.179),
            ('GDP', 'C10H12N5O11P2', -3, 440.178),
            ('AMP', 'C10H12N5O7P', -2, 345.208),
            ('CMP', 'C9H12N3O8P', -2, 321.182),
            ('GMP', 'C10H12N5O8P', -2, 361.207),
            ('UMP', 'C9H11N2O9P', -2, 322.166),
            ('H2O', 'H2O', 0, 18.015),
            ('PPI', 'HO7P2', -3, 174.949),
            ('H', 'H1', 1, 1.008),
            ('P', 'P1', 1, 30.974),
        ]

        for id, formula, charge, molecular_weight in species_types:
            species_type = self.model.species_types.create(id=id, empirical_formula=formula, charge=charge, molecular_weight=molecular_weight)
            species = species_type.species.create(compartment=compartment)
            species.concentration = wc_lang.core.Concentration(value=5e-3, units='M')
