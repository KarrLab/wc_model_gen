""" Generator for protein  degradation submodels based on KBs for random in silico organisms

:Author: Bilal Shaikh <bilal.shaikh@columbia.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
         Jonathan Karr <karr@mssm.edu>
:Date: 2018-07-05
:Copyright: 2018, Karr Lab
:License: MIT
"""
import numpy
import scipy
import wc_kb
import wc_lang
import wc_model_gen
from wc_model_gen.prokaryote.species import SpeciesGenerator


class ProteinDegradationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Gnerator for Protein degradation model"""

    def gen_compartments(self):
        self.cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get_or_create(id='c')
        cytosol.name = 'cytosol'
        cytosol.initial_volume = self.cell.properties.get_one(
            id='initial_volume').value

    def gen_species(self):
        "Generate the protein species for the model"

        speciesGen = SpeciesGenerator(self.knowledge_base, self.model)
        speciesGen.run()

    def gen_reactions(self):

        model = self.model
        submodel = self.submodel
        cytosol = model.compartments.get_one(id='c')

        proteins = self.cell.species_types.get(
            __type=wc_kb.core.ProteinSpeciesType)

        for kb_protein in proteins:
            if kb_protein.id.startswith('protein_'):
                rxn = submodel.reactions.get_or_create(
                    id=kb_protein.id.replace('protein', 'protein_degradation_'))
                rxn.name = kb_protein.name.replace(
                    'protein ', 'protein degradation ')
            else:
                rxn = submodel.reactions.get_or_create(
                    id='protein_degradation_'+str(kb_protein.id))
                rxn.name = 'protein degradation '+str(kb_protein.name)

            model_protein = model.species_types.get_one(
                id=kb_protein.id).species.get_one(compartment=cytosol)

            seq = kb_protein.get_seq()

            rxn.participants = []

            # The protein being degraded
            rxn.participants.add(
                model_protein.species_coefficients.get_or_create(coefficient=-1))

            # ATP used to attach protein to proteosome
            atp = model.species_types.get_one(
                id='atp').species.get_one(compartment=cytosol)
            adp = model.species_types.get_one(
                id='adp').species.get_one(compartment=cytosol)
            pi = model.species_types.get_one(
                id='pi').species.get_one(compartment=cytosol)
            rxn.participants.add(
                atp.species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(
                adp.species_coefficients.get_or_create(coefficient=1))
            rxn.participants.add(
                pi.species_coefficients.get_or_create(coefficient=1))

            # Water needed for the seperation of each amino acid
            h2o = model.species_types.get_one(
                id='h2o').species.get_one(compartment=cytosol)

            rxn.participants.add(
                h2o.species_coefficients.get_or_create(coefficient=-(len(seq)-1)))

            # The 20 amino acids
            amino_acids = ['ala', 'arg', 'asp', 'asn', 'cys', 'gln', 'glu', 'gly', 'his',
                           'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val']
            aas = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                   "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

            for amino_acid, aa in zip(amino_acids, aas):
                species = model.species_types.get_one(
                    id=amino_acid).species.get_one(compartment=cytosol)
                rxn.participants.add(
                    species.species_coefficients.get_or_create(coefficient=seq.count(aa)))

    def gen_rate_laws(self):
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_or_create(id='c')

        proteosome_conc = 5000/scipy.constants.Avogadro / \
            cytosol.initial_volume  # PubMed ID16135238

        deg_protease = model.observables.get_one(id='deg_protease_obs')

        prots = cell.species_types.get(
            __type=wc_kb.ProteinSpeciesType)
        for prot, rxn in zip(prots, self.submodel.reactions):
            rl = rxn.rate_laws.create()
            rl.direction = wc_lang.RateLawDirection.forward

            rl.equation = wc_lang.RateLawEquation(
                expression='{0}[c] * (((k_cat * {1}) / (k_m + {1})) + {2})'.format(prot.id, deg_protease.species[0].species.id(), '0.1'))

            rl.k_cat = 2 * numpy.log(2) / prot.half_life
            rl.k_m = proteosome_conc
            # rl.equation.observables.append(deg_protease)
            rl.equation.modifiers.append(deg_protease.species[0].species)
            rl.equation.modifiers.append(rxn.participants[0].species)
