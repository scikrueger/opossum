#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 15:21:41 2024

@author: flo
"""

import cobra
import pandas as pd
import corpse
import pandas
import networkx as nx


mod = cobra.io.read_sbml_model("/home/flo/Documents/Model_merging/Merged_model_19_08_2023.xml")
mod.solver = "cplex"

# First we get shadow costs before we mess with the models objective
# TODO: Implement handling of custom objective for shadow cost.


raw_diet = pd.read_csv("~/Documents/resources/msb-diets-master/Human/Matjes/matjes-diet30112020.csv")

# reformat diet to fit the model
diet = raw_diet.iloc[1:, [0, 1]]
diet.columns = ['metabolite_id', 'uptake_limit']
# Replace all occurrences of '(e)' with '[f]' in 'metabolite_id'
diet['metabolite_id'] = diet['metabolite_id'].str.replace('(e)', '[f]', regex=False)
diet['uptake_limit'] = diet['uptake_limit'].clip(upper=1000)

# Apply the modified diet constraints to the model
for index, row in diet.iterrows():
    metabolite_id = row['metabolite_id']
    uptake_limit = row['uptake_limit']
    exchange_reaction_id = metabolite_id

    if exchange_reaction_id in mod.reactions:
        mod.reactions.get_by_id(exchange_reaction_id).lower_bound = -uptake_limit
    else:
        print(f"Exchange reaction for {metabolite_id} not found in the model.")
        
#%%        
# =============================================================================
#          DEFINING OBJECTIVE FUNCTION FOR REDUCED COST CALCULATIONS
# =============================================================================
   
# biomass function is defined for reduced cost calculations
# actually unclear what kind of objective function the merged models will have   
# so this is on the backburner right now      

## Approach 1 to set biomass as objective:
biomass_reactions = [reaction for reaction in mod.reactions if 'biomass' in reaction.id.lower() or 'biomass' in reaction.name.lower()]

# Dictionary to store shadow costs for each reaction after each optimization
reduced_costs = pd.Series()

# Loop over each biomass reaction and optimize
for biomass_reaction in biomass_reactions:
    # Set the objective to this biomass reaction
    mod.objective = biomass_reaction
    
    # Optimize the model
    solution = mod.optimize()
    
    # Retrieve shadow costs (dual values)
    reduced_costs = reduced_costs.append(solution.reduced_costs)

## Approach 2 to set biomass as objective:
for reaction in mod.reactions:
    if 'biomass' in reaction.id.lower() or 'biomass' in reaction.name.lower():
        reaction.objective_coefficient = 1
        
mod.objective = sum(reaction.flux_expression * reaction.objective_coefficient for reaction in mod.reactions)
# Idee: Wäre es nicht sinnvoller jeder einzelnt laufen zu lassen und dann summieren?
# Weil so wird vermutlich immer eine einer anderen bevorzugt, was zu einer bevorzugung eines compartments (i.e. 1 bestimmtes Bakterium oder Mikrobiom oder Host führt, je nach dem wie Biomassereaktionen vorhanden und definiert sind)


#%% 
# the model needs nutrients
# ideally it would use existing diet, if none is present use this:
# Problem: How would you identify the custom exchange department?    

# f is the exchange department where nutrients come from    
reactions_in_f = [(reaction.id, reaction.bounds, reaction.compartments) for reaction in mod.reactions if reaction.compartments == {"f"}]




# now get shadow prices
solution = mod.optimize()
shadow_prices = solution.shadow_prices
# now we have them saved and can access them later when deciding which reaction to adapt first
# careful, wee need to consider the absolute, since we want to avoid *any* affect on the objective function
 
# Generate a list of data about reactions confined to the 'f' compartment
data = [{
    'Reaction ID': reaction.id,
    'low_bnd': reaction.lower_bound,
    'upp_bnd': reaction.upper_bound,
} for reaction in mod.reactions if reaction.compartments == {'f'}]

# Create a DataFrame from the list of dictionaries
df = pd.DataFrame(data)



# Access the model's objective
objective_coefficients = {}
for reaction in mod.reactions:
    objective_coefficients[reaction.id] = reaction.objective_coefficient
ob = pandas.DataFrame.from_dict(objective_coefficients, orient='index', columns=['obj_coeff'])


ATP_reactions1 = [[reaction.id, reaction.reaction] for reaction in mod.reactions if "atp" in reaction.reaction]


#%%
# =============================================================================
# # Get all ATP Reactions!
# =============================================================================
# I want all reactions that produce ATP. 
# Problem is that they might be defined in reverse, so just "atp" in products doesnt work,
# and i also want to exclude reactions that just transport ATP, since ATP is not created there
# Ergo: string "atp" is either in .reactants or .products, but NOT in both

# Make it a function for reuse
def getATPreacts():
    # Initialize an empty list to store ATP-producing reactions
    ATP_reactions = []
    
    # Iterate over all reactions in the model
    for reaction in mod.reactions:
        # Check if "atp" appears in the reactants or in the products, but not in both
        reactants_contain_atp = any(str(metabolite).startswith("atp") for metabolite in reaction.reactants)
        products_contain_atp = any(str(metabolite).startswith("atp") for metabolite in reaction.products)
        
        if ((reactants_contain_atp and not products_contain_atp and reaction.lower_bound < 0) or 
            (products_contain_atp and not reactants_contain_atp and reaction.upper_bound > 0)):
            # Append the reaction to the list of ATP-producing reactions
            ATP_reactions.append(reaction)
            
    return(ATP_reactions)


# get the ATP reacts on the virgin model
ATP_reactions = getATPreacts()

# necessary for some reason here
mod.optimize()

# Print the ATP-producing reactions
for reaction in ATP_reactions:
    print(reaction.id, reaction.reaction, reaction.flux)
    
# Get list to feed as objective function
ATPlist = []
for reaction in ATP_reactions:
    ATPlist.append(reaction.id)
    
# update objective functions to ATP producers
mod.objective = ATPlist

# solve mod again with the objective function
mod.optimize()

# update ATP reactions to save new fluxes
ATP_reactions = getATPreacts()



#%%
# close themall exchange reractions off here anyway to make it model independent
# TODO
[[reaction.id, reaction.reaction, reaction.flux] for reaction in mod.reactions if reaction.flux != 0]
# This merged model has custom compartments that need closing

# start with getting all reactions from a compartment (maybe filter only for reaction that involve 2 compartments?)
# now get those
def getCompartmentreacts(comp):
    if (comp in mod.compartments):
        # Use a list comprehension to filter reactions based on compartment
        reactions_in_compartment = [reaction for reaction in mod.reactions 
                                    if any(metabolite.compartment == comp for metabolite in reaction.metabolites)]
        return(reactions_in_compartment)
    
    else: 
        print("The model does not contain such a compartment")
    
f_reacts = getCompartmentreacts("f")
[[reaction.id, reaction.reaction, reaction.flux] for reaction in f_reacts if reaction.flux != 0]
#

# =============================================================================
# Remove all uptake of metabolites
# =============================================================================
# this identifies all exchange reactions independent of compartment
exchange_reactions = [rxn for rxn in mod.reactions if (len(rxn.reactants) == 0 or len(rxn.products) == 0)]

# Create a dictionary to store original bounds
original_bounds = {}

for rxn in exchange_reactions:
    # Store bounds with the reaction ID as the key
    original_bounds[rxn.id] = (rxn.lower_bound, rxn.upper_bound)


# this identifies all exchange reactions independent of compartment
exchange_reactions = [rxn for rxn in mod.reactions if (len(rxn.reactants) == 0 or len(rxn.products) == 0)]

for rxn in exchange_reactions:
    # Check if the reaction is an uptake reaction (no reactants)
    if len(rxn.reactants) == 0:
        rxn.upper_bound = 0  # Block uptake
    # Check if the reaction is a secretion reaction (no products)
    elif len(rxn.products) == 0:
        rxn.lower_bound = 0  # Prevent negative flux (reverse uptake), allow secretion
        
mod.optimize() # mod needs to be solved when bounds are changed

# At the VERY end when all cycles are removed, restore original exchange boundaries
for rxn in exchange_reactions:
    # Retrieve the original bounds from the dictionary
    original_lb, original_ub = original_bounds[rxn.id]
    rxn.lower_bound = original_lb
    rxn.upper_bound = original_ub

#%%
# then close them all off so we can actually see the TICs

ATPflux = 0
for reaction in ATP_reactions:
    ATPflux = ATPflux + reaction.flux
    
________

solution = mod.optimize()

# Threshold for considering flux as non-zero
flux_threshold = 1e-6

# Reactions carrying flux
active_reactions = [
    rxn for rxn in mod.reactions if abs(solution.fluxes[rxn.id]) > flux_threshold
]

print("Reactions involved in ATP-producing cycle:")
for rxn in active_reactions:
    flux = solution.fluxes[rxn.id]
    print(f"{rxn.id}: {rxn.reaction}, Flux: {flux}")

# import networkx as nx

# # Create a directed graph
# G = nx.DiGraph()

# # Add reactions and metabolites as nodes, and connect them
# for rxn in active_reactions:
#     flux = solution.fluxes[rxn.id]
#     for met in rxn.metabolites:
#         coef = rxn.get_coefficient(met)
#         if coef < 0:  # Reactant
#             G.add_edge(met.id, rxn.id)
#         else:  # Product
#             G.add_edge(rxn.id, met.id)

# # Find cycles in the subgraph
# cycles = list(nx.simple_cycles(G))
# print(f"Found {len(cycles)} cycles.")


# # Identify cycles that include ATP metabolite and ATP-producing reactions
# atp_cycles = [cycle for cycle in cycles if any(rxn.id in cycle for rxn in ATP_reactions)] # needs filter for just reactions with atp in it

# print(f"ATP-producing cycles:")
# for i, cycle in enumerate(atp_cycles, 1):
#     print(f"Cycle {i}: {cycle}")


#the cycles above are just 1 metabolite 
#ALTERNATIVE APPROACH, reactions only###

import networkx as nx

# Assume 'model' is your cobra model object

# Create a directed graph of reactions only
G_reactions = nx.DiGraph()

# Add all reactions as nodes
for rxn in model.reactions:
    G_reactions.add_node(rxn.id)

# Create edges: R_i -> R_j if R_i produces something R_j consumes
for rxn_i in model.reactions:
    # Get the set of products from rxn_i
    products_i = {met.id for met in rxn_i.products}

    for rxn_j in model.reactions:
        if rxn_i.id == rxn_j.id:
            # Skip self-loops for clarity (can remove if you want to detect them)
            continue
        
        # Get the reactants of rxn_j
        reactants_j = {met.id for met in rxn_j.reactants}

        # If there's an intersection, it means rxn_i leads to a metabolite rxn_j needs
        if products_i & reactants_j:
            # Add an edge from rxn_i to rxn_j
            G_reactions.add_edge(rxn_i.id, rxn_j.id)

# Now detect cycles in the reaction-only graph
reaction_cycles = list(nx.simple_cycles(G_reactions))

print(f"Found {len(reaction_cycles)} reaction-level cycles.")
for i, cyc in enumerate(reaction_cycles, start=1):
    print(f"Cycle {i}: {cyc}")


from collections import defaultdict

for i, cyc in enumerate(reaction_cycles, start=1):
    # Stoichiometry balance
    met_balance = defaultdict(float)

    # Aggregate the stoichiometry of all reactions in the cycle
    for rxn_id in cyc:
        rxn = model.reactions.get_by_id(rxn_id)
        for met, coeff in rxn.metabolites.items():
            met_balance[met.id] += coeff

    # Check if all metabolites (except perhaps ATP/ADP) balance out
    # A true futile cycle would have zero net consumption for all non-energy metabolites
    # and positive net production of ATP.
    all_balanced = True
    for met_id, net_coeff in met_balance.items():
        # If there's a net consumption of a metabolite that can't be produced internally, it's not a closed cycle.
        # If the cycle is supposed to be infeasible, you’ll see that something doesn’t balance out.
        # Decide on rules: for a pure "futile cycle", you'd expect net zero for non-energy metabolites.
        if met_id not in ['atp_c', 'adp_c', 'pi_c'] and abs(net_coeff) > 1e-6:
            all_balanced = False
            break

    if all_balanced:
        print(f"Cycle {i} is balanced and could represent a futile cycle: {cyc}")
    else:
        print(f"Cycle {i} is not fully balanced: {cyc}")



#TODO:
    # build proper nodes with only fully accounted reactions for all reactions that have flux.
    
    
    import networkx as nx

G_reactions = nx.DiGraph()

for rxn in model.reactions:
    G_reactions.add_node(rxn.id)

for rxn_i in model.reactions:
    # Extract products of rxn_i
    products_i = {met.id: coeff for met, coeff in rxn_i.metabolites.items() if coeff > 0}
    
    for rxn_j in model.reactions:
        if rxn_i == rxn_j:
            continue
        
        # Extract reactants of rxn_j (coeff < 0)
        reactants_j = {met.id: -coeff for met, coeff in rxn_j.metabolites.items() if coeff < 0}

        # Check if rxn_i supplies all reactants_j in necessary amounts
        can_supply_all = True
        for met_id, needed_amount in reactants_j.items():
            produced_amount = products_i.get(met_id, 0)
            if produced_amount < needed_amount:
                can_supply_all = False
                break

        # Add edge only if full stoichiometric requirements are met
        if can_supply_all:
            G_reactions.add_edge(rxn_i.id, rxn_j.id)

# Now detect cycles
reaction_cycles = list(nx.simple_cycles(G_reactions))

print(f"Found {len(reaction_cycles)} reaction-level cycles with full substrate coverage:")
for i, cyc in enumerate(reaction_cycles, start=1):
    print(f"Cycle {i}: {cyc}")

###########################################################################################
###################################   eQuilibrium   #######################################
###########################################################################################
from equilibrator_api import Reaction, ComponentContribution, Q_, default_compound_matcher
from equilibrator_api import CompoundMatcher
cc = ComponentContribution()

atpase_reaction = cc.parse_reaction_formula("chebi:CHEBI:30616 + kegg:C00001 = hmdb:HMDB01341 + seed:cpd00009")

standard_dg_prime = cc.standard_dg_prime(atpase_reaction)
physiological_dg_prime = cc.physiological_dg_prime(atpase_reaction)


test1 = reaction_to_bigg_string(active_reactions[0])
test1 = 'h2o[c] + datp[c] + trdox[c] = atp[c] + trdrd[c]' #doesnt work
test1 = 'bigg:h2o + bigg:datp + bigg:trdox = bigg:atp + bigg:trdrd' #doesnt work
test1 = 'bigg.metabolite:h2o[c] + bigg.metabolite:datp[c] + bigg.metabolite:trdox[c] = bigg.metabolite:atp[c] + bigg.metabolite:trdrd[c]' #doesnt work
test1 = 'bigg.metabolite:h2o + bigg.metabolite:datp + bigg.metabolite:trdox = bigg.metabolite:atp + bigg.metabolite:trdrd' #works
test1 = 'h2o + datp + trdox = atp + trdrd' #works
test2 = cc.parse_reaction_formula(test1)
cc.standard_dg_prime(test2)

def reaction_to_bigg_string(rxn):
    left = []
    right = []
    for met, coeff in rxn.metabolites.items():
        name = f"bigg.metabolite:{met.id}"
        if coeff < 0:
            left.append(f"{-coeff:g} {name}" if abs(coeff) != 1 else f"{name}")
        elif coeff > 0:
            right.append(f"{coeff:g} {name}" if abs(coeff) != 1 else f"{name}")
    return ' + '.join(left) + ' = ' + ' + '.join(right)

for reaction in active_reactions:
    rxn_str = reaction_to_bigg_string(reaction)
    try:
        # Pass both arguments explicitly!
        eq_rxn = cc.parse_reaction_formula(rxn_str)
        dG0 = cc.standard_dg_prime(eq_rxn)
        print(f"{reaction.id}: ΔG'° = {dG0.value:.2f} ± {dG0.uncertainty:.2f} {dG0.units}")
    except Exception as e:
        print(f"{reaction.id}: ΔG calculation failed: {e}")



# free energ calc:
# cycle through reactions
# access their compounds with .products and .inputs (other name)

# for now treat everything as single compartment (or exclude multicompartment cases for now and use simpler approach??)
# eQuilibrator supports multicompartment, but needs more bio-data for that
# reformat them and calc deltaG with eQuilibrium
# gather results in dataframe

# pFBA for shadow costs (potentially right at the beginning before objective function is changed to ATP)


