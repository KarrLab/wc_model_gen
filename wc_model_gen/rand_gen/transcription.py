""" Generator for transcription submodels based on KBs for random in silico organisms

:Author: Jonathan Karr <karr@mssm.edu>
         Ashwin Srinivasan
:Date: 2018-06-11
:Copyright: 2018, Karr Lab
:License: MIT
"""

import numpy
import scipy
import wc_kb
import wc_lang
import wc_model_gen


class TranscriptionSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for transcription submodel """

    def gen_compartments(self):
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get_or_create(id='c')
        cytosol.name = 'cytosol'
        cytosol.initial_volume = cell.properties.get_one(id='mean_volume').value

    def gen_species(self):
        """ Generate species associated with submodel """
        cell = self.knowledge_base.cell
        model = self.model
        cytosol = model.compartments.get(id='c')[0]

         # get or create RNA species
        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna in rnas:
            species_type = model.species_types.get_or_create(id=rna.id)
            if not species_type.name:
                species_type.name = rna.name
                species_type.type = wc_lang.SpeciesTypeType.rna
                species_type.structure = rna.get_seq()
                species_type.empirical_formula = rna.get_empirical_formula()
                species_type.molecular_weight = rna.get_mol_wt()
                species_type.charge = rna.get_charge()
                species = species_type.species.get_or_create(
                    compartment=cytosol)
                species.concentration = wc_lang.Concentration(
                    value=rna.concentration, units=wc_lang.ConcentrationUnit.M)

        for protein in self.knowledge_base.cell.species_types.get(__type=wc_kb.core.ProteinSpeciesType):

            species_type = self.model.species_types.get_or_create(
                id=protein.id)
            if not species_type.name:
                # Add functional form of protein
                species_type.name = protein.name
                species_type.type = wc_lang.SpeciesTypeType.protein
                species_type.structure = protein.get_seq()
                species_type.empirical_formula = protein.get_empirical_formula()
                species_type.molecular_weight = protein.get_mol_wt()
                species_type.charge = protein.get_charge()
                species = species_type.species.get_or_create(
                    compartment=cytosol)

                species.concentration = wc_lang.Concentration(
                    value=protein.concentration, units=wc_lang.ConcentrationUnit.M)

        for kb_observable in self.knowledge_base.cell.observables:
            model_observable = self.model.observables.get_or_create(
                id=kb_observable.id)
            if not model_observable.name:
                model_observable.name = kb_observable.name
                for kb_species_coefficient in kb_observable.species:
                    kb_species = kb_species_coefficient.species
                    kb_species_type = kb_species.species_type
                    kb_compartment = kb_species.compartment
                    model_species_type = model.species_types.get_one(
                        id=kb_species_type.id)
                    model_species = model_species_type.species.get_one(
                        compartment=model.compartments.get_one(id=kb_compartment.id))
                    model_coefficient = kb_species_coefficient.coefficient
                    model_species_coefficient = wc_lang.SpeciesCoefficient()
                    model_species_coefficient.species = model_species
                    model_species_coefficient.coefficient = model_coefficient

                    model_observable.species.append(model_species_coefficient)
                    
                    model_species.species_coefficients.get_or_create(coefficient = model_coefficient)

                for kb_observable_observable in kb_observable.observables:
                    model_observable_observable = model.observables.get_or_create(
                        id=kb_observable_observable.id)
                    model_observable.observables.append(
                        model_observable_observable)

                
    def gen_reactions(self):
        """ Generate reactions associated with submodel """
        model = self.model
        submodel = self.submodel
        cytosol = model.compartments.get_one(id='c')
        atp = model.species_types.get_one(id='atp').species.get_one(compartment=cytosol)
        ctp = model.species_types.get_one(id='ctp').species.get_one(compartment=cytosol)
        gtp = model.species_types.get_one(id='gtp').species.get_one(compartment=cytosol)
        utp = model.species_types.get_one(id='utp').species.get_one(compartment=cytosol)
        ppi = model.species_types.get_one(id='ppi').species.get_one(compartment=cytosol)
        h2o = model.species_types.get_one(id='h2o').species.get_one(compartment=cytosol)
        h = model.species_types.get_one(id='h').species.get_one(compartment=cytosol)

        cell = self.knowledge_base.cell
        kb_rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for kb_rna in kb_rnas:
            if kb_rna.id.startswith('rna_'):
                rxn = submodel.reactions.get_or_create(id=kb_rna.id.replace('rna_', 'transcription_'))
                rxn.name = kb_rna.name.replace('RNA ', 'Transcription')
            else:
                rxn = submodel.reactions.get_or_create(id='transcription_'+str(kb_rna.id))
                rxn.name = 'Transcription '+str(kb_rna.name)


            model_rna = model.species_types.get_one(id=kb_rna.id).species.get_one(compartment=cytosol)
            seq = kb_rna.get_seq()
            rxn.participants = []
            rxn.participants.add(atp.species_coefficients.get_or_create(coefficient=-seq.count('A')))
            rxn.participants.add(ctp.species_coefficients.get_or_create(coefficient=-seq.count('C')))
            rxn.participants.add(gtp.species_coefficients.get_or_create(coefficient=-seq.count('G')))
            rxn.participants.add(utp.species_coefficients.get_or_create(coefficient=-seq.count('U')))
            rxn.participants.add(h.species_coefficients.get_or_create(coefficient=-(kb_rna.get_len() - 1)))
            rxn.participants.add(model_rna.species_coefficients.get_or_create(coefficient=1))
            rxn.participants.add(ppi.species_coefficients.get_or_create(coefficient=kb_rna.get_len()))
            rxn.participants.add(h2o.species_coefficients.get_or_create(coefficient=kb_rna.get_len() - 1))

    def gen_rate_laws(self):
        """ Generate rate laws associated with submodel """
        model = self.model
        cell = self.knowledge_base.cell
        cytosol = model.compartments.get_one(id='c')

        mean_volume = cell.properties.get_one(id='mean_volume').value
        mean_doubling_time = cell.properties.get_one(id='mean_doubling_time').value
        poly_avg_conc = 3000/scipy.constants.Avogadro / cytosol.initial_volume #http://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=2&id=106199
        rna_poly = self.model.observables.get_one(
            id='rna_poly_obs')
        exp='(((k_cat * {}) / (k_m + {})))'.format(rna_poly.id, rna_poly.id)
        equation = wc_lang.RateLawEquation(expression = exp)
        equation.parameters.append(
                rna_poly)
        
        rnas = cell.species_types.get(__type=wc_kb.RnaSpeciesType)
        for rna, rxn in zip(rnas, self.submodel.reactions):
            rl = rxn.rate_laws.create()
            rl.equation = equation
            rl.direction = wc_lang.RateLawDirection.forward
            rl.k_cat = 2 * (numpy.log(2) / rna.half_life + numpy.log(2) / mean_doubling_time)
            rl.k_m = poly_avg_conc
