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
            rxn = submodel.reactions.get_or_create(id=protein_kb.id.replace('prot_gene_', 'degrad_prot_'))
            rxn.name = protein_kb.id.replace('prot_gene_', 'degrad_prot_')
            rxn.participants = []

            # Adding participants to LHS
            rxn.participants.add(protein_model.species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(atp.species_coefficients.get_or_create(coefficient=-1))
            rxn.participants.add(h2o.species_coefficients.get_or_create(coefficient=-(len(seq)-1)))

            # Adding participants to RHS
            rxn.participants.add(adp.species_coefficients.get_or_create(coefficient=1))
            rxn.participants.add(pi.species_coefficients.get_or_create(coefficient=1))

            # The code below should be used as currently tRNAs and AAs are always associated
            codons=[]
            for start_position in range(0,len(protein_kb.gene.get_seq())-3,3):
                codons.append(str(protein_kb.gene.get_seq()[start_position:start_position+3]))

            for codon in set(codons):
                obs_model_id = 'tRNA_'+ codon + '_obs'
                obs_model = model.observables.get_one(id=obs_model_id)
                for specie in obs_model.expression.species:
                    rxn.participants.add(
                        specie.species_coefficients.get_or_create(coefficient=codons.count(codon)))

            #for amino_acid, aa in zip(amino_acids, aas):
            #    species = model.species_types.get_one(id=amino_acid).species.get_one(compartment=cytosol)
            #    rxn.participants.add(species.species_coefficients.get_or_create(coefficient=seq.count(aa)))

            # Add members of the degradosome
            # Counterintuitively .specie is a KB species_coefficient object
            for degradosome_kb in cell.observables.get_one(id='degrade_protease_obs').species:
                degradosome_species_type_model = model.species_types.get_one(id=degradosome_kb.species.species_type.id)
                degradosome_species_model = degradosome_species_type_model.species.get_one(compartment=cytosol)

                rxn.participants.add(degradosome_species_model.species_coefficients.get_or_create(coefficient=(-1)*degradosome_kb.coefficient))
                rxn.participants.add(degradosome_species_model.species_coefficients.get_or_create(coefficient=degradosome_kb.coefficient))

    def gen_phenomenological_rates(self):
        """ Generate rate laws with exponential dynamics """

        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='protein_degradation')

        proteins_kb = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb, rxn in zip(proteins_kb, submodel.reactions):
            protein_model = model.species_types.get_one(id=protein_kb.id).species[0]
            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = '({} / {}) * {}'.format(numpy.log(2), protein_kb.half_life, protein_model.id())
            rate_law.equation = wc_lang.RateLawEquation(expression = expression)
            rate_law.equation.modifiers.append(rxn.participants[0].species)

    def gen_mechanistic_rates(self):
        """ Generate rate laws associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        submodel = model.submodels.get_one(id='protein_degradation')
        mean_volume = cell.properties.get_one(id='initial_volume').value
        mean_cell_cycle_length = cell.properties.get_one(id='cell_cycle_length').value
        cytosol = cell.compartments.get_one(id='c')

        proteins_kbs = cell.species_types.get(__type=wc_kb.prokaryote_schema.ProteinSpeciesType)
        for protein_kb, rxn in zip(proteins_kbs, submodel.reactions):

            protein_model = model.species_types.get_one(id=protein_kb.id).species[0]
            rate_law = rxn.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            expression = 'k_cat*'
            modifiers = []
            rate_avg = ''
            beta = 2

            #TODO: replace with calculation of avg half life; 553s is avg of Mycoplasma RNAs
            if protein_kb.half_life == 0:
                protein_kb.half_life = 12*60*60

            for participant in rxn.participants:
                if participant.coefficient < 0:
                    avg_conc = (3/2)*participant.species.concentration.value
                    modifiers.append(participant.species)
                    rate_avg += '({}/({}+({}*{})))*'.format(avg_conc, avg_conc, beta, avg_conc)
                    expression += '({}/({}+(3/2)*{}))*'.format(participant.species.id(),
                                                              participant.species.id(),
                                                              participant.species.concentration.value)

            # Clip off trailing * character
            expression = expression[:-1]
            rate_avg = rate_avg[:-1]

            # Create / add rate law equation
            if 'rate_law_equation' not in locals():
                rate_law_equation = wc_lang.RateLawEquation(expression=expression, modifiers=modifiers)

            rate_law.equation = rate_law_equation

            # Calculate k_cat
            exp_expression = '({}*(1/{}+1/{})*{})'.format(
                                numpy.log(2),
                                cell.properties.get_one(id='cell_cycle_length').value,
                                protein_kb.half_life,
                                3/2*protein_kb.species.get_one(compartment=cytosol).concentrations.value)
                                #This should have units of M

            rate_law.k_cat = eval(exp_expression) / eval(rate_avg)
