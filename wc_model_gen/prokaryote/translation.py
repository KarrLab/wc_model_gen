""" Generating wc_lang formatted models from knowledge base.
:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
         Ashwin Srinivasan <ashwins@mit.edu>
:Date: 2018-01-21
:Copyright: 2018, Karr Lab
:License: MIT

TODO: aminoacyl-tRnas are missign from observables
"""

import wc_kb
import wc_lang
import wc_model_gen
import numpy
import scipy
from wc_model_gen.prokaryote.species import SpeciesGenerator


class TranslationSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generate translation submodel. """

    def gen_species(self):
        """ Generate protein species """
        speciesGen = SpeciesGenerator(self.knowledge_base, self.model)
        speciesGen.run()

    def gen_reactions(self):
        """ Generate a lumped reaction that cvers initiation, elongation and termination for each protein translated """
        submodel = self.submodel
        compartment = self.model.compartments.get_one(id='c')

        proteins = self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType)

        # Add translation initating reactions
        for protein_kb in proteins:
            reaction = submodel.reactions.get_or_create(id='translation_' + protein_kb.id)
            reaction.name = 'translate '+ protein_kb.id
            n_steps = protein_kb.get_len()

            # Adding reaction participants LHS
            ribosome = self.model.observables.get_one(id='complex_70S_obs').expression.species[0]
            gtp = self.model.species_types.get_one(id='gtp').species.get_one(compartment=compartment)
            initiation_factors = self.model.observables.get_one(id='translation_init_factors_obs').expression.species[0]
            elongation_factors = self.model.observables.get_one(id='translation_elongation_factors_obs').expression.species[0]
            release_factors = self.model.observables.get_one(id='translation_release_factors_obs').expression.species[0]

            reaction.participants.add(ribosome.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.add(gtp.species_coefficients.get_or_create(coefficient=-(n_steps+2)))
            reaction.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=-1))
            reaction.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=-n_steps))
            reaction.participants.add(release_factors.species_coefficients.get_or_create(coefficient=-1))

            # Add tRNAs to
            bases = "TCAG"
            codons = [a + b + c for a in bases for b in bases for c in bases]

            for codon in codons:
                if codon not in ['TAG', 'TAA', 'TGA']:
                    n = str(protein_kb.gene.get_seq()).count(codon)
                    if n > 0:
                        trna = self.model.observables.get_one(id='tRNA_'+codon+'_obs').expression.species[0]
                        reaction.participants.add(trna.species_coefficients.get_or_create(coefficient=-n))

            # Adding reaction participants RHS
            protein = self.model.species_types.get_one(id=protein_kb.id).species.get_one(compartment=compartment)
            gdp = self.model.species_types.get_one(id='gdp').species.get_one(compartment=compartment)
            pi = self.model.species_types.get_one(id='pi').species.get_one(compartment=compartment)

            #ribosome = self.model.observables.get_one(id='complex_70S_obs').expression.species[0]
            #initiation_factors = self.model.observables.get_one(id='translation_init_factors_obs').expression.species[0]
            #elongation_factors = self.model.observables.get_one(id='translation_elongation_factors_obs').expression.species[0]
            #release_factors = self.model.observables.get_one(id='translation_release_factors_obs').expression.species[0]

            reaction.participants.add(protein.species_coefficients.get_or_create(coefficient=1))
            reaction.participants.add(ribosome.species_coefficients.get_or_create(coefficient=1))
            reaction.participants.add(gdp.species_coefficients.get_or_create(coefficient=n_steps+2))
            reaction.participants.add(pi.species_coefficients.get_or_create(coefficient=2*n_steps))
            reaction.participants.add(initiation_factors.species_coefficients.get_or_create(coefficient=1))
            reaction.participants.add(elongation_factors.species_coefficients.get_or_create(coefficient=n_steps))
            reaction.participants.add(release_factors.species_coefficients.get_or_create(coefficient=1))

    def gen_rate_laws(self):
        """ Generate rate laws associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='translation')

        mean_volume = cell.properties.get_one(id='initial_volume').value
        mean_doubling_time = cell.properties.get_one(id='doubling_time').value

        # check ribosome count
        poly_avg_conc = 5000/scipy.constants.Avogadro / cytosol.initial_volume
        rna_polymerase = model.observables.get_one(id='rna_polymerase_obs')

        for reaction in submodel.reactions:
            rate_law = reaction.rate_laws.create()
            rate_law.direction = wc_lang.RateLawDirection.forward
            rate_law_equation = wc_lang.RateLawEquation()
            expression = 'k_cat*'

            for participant in reaction.participants:
                if participant.coefficient < 0:
                    rate_law_equation.modifiers.append(participant.species)
                    expression += '({}/({}+(3/2)*{}))*'.format(participant.species.id(),
                                                              participant.species.id(),
                                                              participant.species.concentration.value)

            expression = expression[:-1] #clip off trailing * character
            rate_law_equation.expression = expression
            rate_law.equation = rate_law_equation

        """ Generate the rate laws associate dwith reactions
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')
        submodel = model.submodels.get_one(id='translation')

        mean_volume = cell.properties.get_one(id='initial_volume').value
        mean_doubling_time = cell.properties.get_one(id='doubling_time').value

        IF = self.model.observables.get_one(id='translation_init_factors_obs')
        IF_avg_conc = 0.05  # placeholder
        EF = self.model.observables.get_one(id='translation_elongation_factors_obs')
        EF_avg_conc = 0.05  # placeholder
        RF = self.model.observables.get_one(id='translation_release_factors_obs')
        RF_avg_conc = 0.05  # placeholder

        exp = 'k_cat'

        init_eq = wc_lang.core.RateLawEquation(expression=exp + ' * ({}'.format(IF.expression.species[0].id()) +
                                               '/ (k_m +{}))'.format(IF.expression.species[0].id()))
        init_eq.modifiers.append(IF.expression.species[0])

        elon_eq = wc_lang.core.RateLawEquation(expression=exp + ' * ({}'.format(EF.expression.species[0].id()) +
                                               '/ (k_m +{}))'.format(EF.expression.species[0].id()))
        # elon_eq.observables.append(EF)
        elon_eq.modifiers.append(EF.expression.species[0])

        term_eq = wc_lang.core.RateLawEquation(expression=exp + ' * ({}'.format(RF.expression.species[0].id()) +
                                               '/ (k_m +{}))'.format(RF.expression.species[0].id()))
        # term_eq.observables.append(RF)
        term_eq.modifiers.append(RF.expression.species[0])

        for reaction in self.submodel.reactions:
            if reaction.id.startswith('translation_init_'):
                rl = reaction.rate_laws.create()
                rl.direction = wc_lang.RateLawDirection.forward
                prot_id = reaction.id[reaction.id.find('init_')+5:]
                prot = cell.species_types.get_one(id=prot_id)
                rl.k_cat = 2 * (numpy.log(2) / prot.half_life +
                                numpy.log(2) / mean_doubling_time)
                rl.equation = init_eq
                rl.k_m = IF_avg_conc
            elif reaction.id.startswith('translation_elon_'):
                rl = reaction.rate_laws.create()
                rl.direction = wc_lang.RateLawDirection.forward
                prot_id = reaction.id[reaction.id.find('elon_')+5:]
                prot = cell.species_types.get_one(id=prot_id)
                rl.k_cat = 2 * (numpy.log(2) / prot.half_life +
                                numpy.log(2) / mean_doubling_time)
                rl.equation = elon_eq
                rl.k_m = EF_avg_conc

            else:
                rl = reaction.rate_laws.create()
                rl.direction = wc_lang.RateLawDirection.forward
                prot_id = reaction.id[reaction.id.find('term_')+5:]
                prot = cell.species_types.get_one(id=prot_id)
                rl.k_cat = 2 * (numpy.log(2) / prot.half_life +
                                numpy.log(2) / mean_doubling_time)
                rl.equation = term_eq
                rl.k_m = RF_avg_conc
                """
