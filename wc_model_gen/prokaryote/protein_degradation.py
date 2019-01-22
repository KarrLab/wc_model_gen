""" Generator for protein  degradation submodels based on KBs for random in silico organisms

:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
         Jonathan Karr <karr@mssm.edu>
:Date: 2018-07-05
:Copyright: 2018, Karr Lab
:License: MIT
"""

from wc_utils.util.units import unit_registry
import wc_model_gen
import wc_lang
import wc_kb
import numpy
import math


class ProteinDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for protein degradation model"""

    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='protein_degradation')
        cytosol = model.compartments.get_one(id='c')

        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        adp = model.species_types.get_one(id='adp').species.get_one(compartment=cytosol)
        pi = model.species_types.get_one(id='pi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)

        amino_acids = ['ala', 'arg', 'asp', 'asn', 'cys', 'gln', 'glu', 'gly', 'his',
                       'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val']

        aas = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
               "S", "T", "W", "Y", "V"]

        proteins_kb = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)

        for protein_kb in proteins_kb:

            protein_model = model.species_types.get_one(id=protein_kb.id).species.get_one(compartment=cytosol)
            seq = protein_kb.get_seq()
            reaction = model.reactions.get_or_create(submodel=submodel, id='degradation_' + protein_kb.id)
            reaction.name = 'degradation ' + protein_kb.name
            reaction.participants = []

            # Adding participants to LHS
            reaction.participants.add(protein_model.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.add(atp.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.add(h2o.species_coefficients.get_or_create(coefficient=-(len(seq)-1)))

            # Adding participants to RHS
            reaction.participants.add(adp.species_coefficients.get_or_create(coefficient=1))
            reaction.participants.add(pi.species_coefficients.get_or_create(coefficient=1))

            # The code below should be used as currently tRNAs and AAs are always associated
            codons = []
            for start_position in range(0, len(protein_kb.gene.get_seq())-3, 3):
                codons.append(str(protein_kb.gene.get_seq()[start_position:start_position+3]))

            for codon in set(codons):
                obs_model_id = 'tRNA_' + codon + '_obs'
                obs_model = model.observables.get_one(id=obs_model_id)
                for specie in obs_model.expression.species:
                    reaction.participants.add(
                        specie.species_coefficients.get_or_create(coefficient=codons.count(codon)))

            # for amino_acid, aa in zip(amino_acids, aas):
            #    species = model.species_types.get_one(id=amino_acid).species.get_one(compartment=cytosol)
            #    reaction.participants.add(species.species_coefficients.get_or_create(coefficient=seq.count(aa)))

            # Add members of the degradosome
            # Counterintuitively .specie is a KB species_coefficient object
            for degradosome_kb in cell.observables.get_one(id='degrade_protease_obs').expression.species:
                degradosome_species_type_model = model.species_types.get_one(id=degradosome_kb.species_type.id)
                degradosome_species_model = degradosome_species_type_model.species.get_one(compartment=cytosol)

                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(
                    coefficient=-1))
                reaction.participants.add(degradosome_species_model.species_coefficients.get_or_create(
                    coefficient=1))

    def gen_phenom_rates(self):
        """ Generate rate laws with exponential dynamics """
        model = self.model
        kb = self.knowledge_base
        submodel = model.submodels.get_one(id='protein_degradation')
        cytosol = model.compartments.get_one(id='c')
        proteins_kb = kb.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        mean_doubling_time = kb.cell.properties.get_one(id='mean_doubling_time').value

        for protein_kb, reaction in zip(proteins_kb, self.submodel.reactions):
            species_type_model = model.species_types.get_one(id=protein_kb.id)
            species_model = species_type_model.species.get_one(compartment=cytosol)

            rate_law = model.rate_laws.create(
                reaction=reaction,
                direction=wc_lang.RateLawDirection.forward)
            rate_law.id = rate_law.gen_id()

            half_life_model = model.parameters.get_or_create(id='half_life_{}'.format(protein_kb.id),
                                                             type=None,
                                                             value=protein_kb.half_life,
                                                             units=unit_registry.parse_units('s'))
            molecule_units = model.parameters.get_or_create(id='molecule_units',
                                                            type=None,
                                                            value=1.,
                                                            units=unit_registry.parse_units('molecule'))

            expression = '(log(2) / {}) / {} * {}'.format(half_life_model.id, molecule_units.id, species_model.id)

            objects = {
                wc_lang.Parameter: {half_life_model.id: half_life_model, molecule_units.id: molecule_units},
                wc_lang.Species: {species_model.id: species_model},
            }
            rate_law.expression, error = wc_lang.RateLawExpression.deserialize(expression, objects)
            assert error is None, str(error)

    def gen_mechanistic_rates(self):
        """ Generate rate laws associated with submodel """
        submodel = self.model.submodels.get_one(id='protein_degradation')
        proteins_kb = self.knowledge_base.cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        mean_doubling_time = self.knowledge_base.cell.properties.get_one(id='mean_doubling_time').value

        for protein_kb, reaction in zip(proteins_kb, self.submodel.reactions):
            self.gen_mechanistic_rate_law_eq(specie_type_kb=protein_kb,
                                             submodel=submodel,
                                             reaction=reaction,
                                             beta=1,
                                             half_life=protein_kb.half_life,
                                             mean_doubling_time=mean_doubling_time)
