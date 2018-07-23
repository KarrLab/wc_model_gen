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

            # Add amino acids
            # CHARGE VALUES ARE PLACEHOLDERS! - IT IS PH DEPENDENT, NEED TO RESERACH BIOCHEM
            ('Ala', 'C3H7NO2', 0, 89.0935), #Alanine
            ('Arg', 'C6H14N4O2', 0, 174.2017), #Arginine
            ('Asn', 'C4H8N2O3', 0, 132.1184), #Asparagine
            ('Asp', 'C4H7NO4', 0, 133.1032), #Aspartic acid
            ('Cys', 'C3H7NO2S', 0, 121.1590), #Cysteine
            ('Gln', 'C5H10N2O3', 0, 147.1299), #Glutamine
            ('Glu', 'C5H9NO4',  0, 146.1451), #Glutamic acid
            ('Gly', 'C2H5NO2', 0, 75.0669), #Glycine
            ('His', 'C6H9N3O2', 0, 155.1552), #Histidine
            ('Ile', 'C6H13NO2', 0, 131.1736), #Isoleucine
            ('Leu', 'C6H13NO2', 0, 131.1736), #Leucine
            ('Lys', 'C6H14N2O2', 0, 146.1882), #Lysine
            ('Met', 'C5H11NO2S',  0, 149.2124), #Methionine
            ('Phe', 'C9H11NO2', 0, 165.1900), #Phenylalanine
            ('Pro', 'C5H9NO2', 0, 115.1310), #Proline
            ('Ser', 'C3H7NO3', 0, 105.0930), #Serine
            ('Thr', 'C4H9NO3', 0, 119.1197), #Threonine
            ('Trp', 'C11H12N2O2', 0, 204.2262), #Tryptophan
            ('Tyr', 'C9H11NO3', 0, 181.1894), #Tyrosine
            ('Val', 'C5H11NO2', 0, 117.1469), #Valine

            ('H2O', 'H2O', 0, 18.015),
            ('PPI', 'HO7P2', -3, 174.949),
            ('H', 'H1', 1, 1.008),
            ('P', 'P1', 1, 30.974),
            ('N', 'N1', 1, 14.0067),
            ('C', 'C1', 1, 12.0107),
            ('O', 'O1', 1, 15.994)
        ]

        for id, formula, charge, molecular_weight in species_types:
            species_type = self.model.species_types.create(id=id, empirical_formula=formula, charge=charge, molecular_weight=molecular_weight)
            specie = species_type.species.create(compartment=compartment)
            specie.concentration = wc_lang.core.Concentration(value=0.005, units=wc_lang.ConcentrationUnit.uM)

        # Add water to extracellular space so compartment 'e' does not have 0 mass
        compartment = self.model.compartments.get_one(id='e')
        specie = self.model.species_types.get_one(id='H2O').species.create(compartment=compartment)
        specie.concentration = wc_lang.core.Concentration(value=2, units=wc_lang.ConcentrationUnit.uM)

        # Add extra water to intracellular space
        compartment = self.model.compartments.get_one(id='c')
        specie = self.model.species_types.get_one(id='H2O').species.get_one(compartment=compartment)
        specie.concentration.value = 0.05
