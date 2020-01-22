""" Generator for metabolism submodel for eukaryotes

:Author: Yin Hoon Chew <yinhoon.chew@mssm.edu>
:Date: 2020-01-21
:Copyright: 2020, Karr Lab
:License: MIT
"""

import wc_model_gen


class MetabolismSubmodelGenerator(wc_model_gen.SubmodelGenerator):
    """ Generator for metabolism submodel 

        Options:
        * recycled_metabolites (:obj:`dict`): a dictionary with metabolite species IDs
            as keys and recycled amounts as values
        * atp_from_submodel (:obj:`bool`): if True, ATP requirement will be calculated 
            from other generated submodels, else the value will be expected to be in a 
            parameter created during model initialization; default is True     
    """

    def clean_and_validate_options(self):
        """ Apply default options and validate options """
        options = self.options

        recycled_metabolites = options.get('recycled_metabolites', {})
        options['recycled_metabolites'] = recycled_metabolites

        atp_from_submodel = options.get('atp_from_submodel', True)
        options['atp_from_submodel'] = atp_from_submodel

    def gen_reactions(self):
        """ Generate reactions associated with submodel 
        
        A biomass reaction is generated by accounting for all the metabolites that 
        are consumed and produced by the reactions in other submodels, and the 
        metabolites that are in the free pool

        """
        cell = self.knowledge_base.cell     
        model = self.model
        submodel = self.submodel
        
        biomass_rxn = model.reactions.create(
            submodel=submodel,
            id='biomass_reaction',
            name='Biomass reaction',
            reversible=False,
            comments='Pseudo-reaction for use as dFBA objective function')

        # Add metabolites in the free pool to the RHS of the biomass reaction
        biomass_rxn.participants.add(
            model_species.species_coefficients.get_or_create(coefficient=participant.coefficient))

        # Add metabolites to be recycled to the LHS of the biomass reaction
        
        # Add metabolite requirement of other submodels
		# Add ATP usage

        # Add biomass reaction as objective function
        submodel.dfba_obj = wc_lang.DfbaObjective(model=model)
        submodel.dfba_obj.id = submodel.dfba_obj.gen_id()
        obj_expression = biomass_rxn.id
        dfba_obj_expression, error = wc_lang.DfbaObjectiveExpression.deserialize(
            obj_expression, {wc_lang.Reaction: {biomass_rxn.id: biomass_rxn}})
        assert error is None, str(error)
        submodel.dfba_obj.expression = dfba_obj_expression

    def gen_rate_laws(self):
        """ Generate rate laws for the reactions in the submodel """
        pass
        
    def calibrate_submodel(self):
        """ Calibrate the submodel using data in the KB """
        pass