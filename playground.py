#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 16:01:50 2024

@author: flo
"""

import cobra
from cobra import Model, Reaction, Metabolite

# Create a new model
model = Model('test_model')

# Create metabolites
A = Metabolite('A', compartment='c')
B = Metabolite('B', compartment='c')
C = Metabolite('C', compartment='c')
D = Metabolite('D', compartment='c')
E = Metabolite('E', compartment='c')
ATP = Metabolite('atp_c', compartment='c')
ADP = Metabolite('adp_c', compartment='c')
P_i = Metabolite('pi_c', compartment='c')

# Create reactions
R1 = Reaction('R1')
R1.name = 'Reaction 1'
R1.lower_bound = 0
R1.upper_bound = 1000
R1.add_metabolites({A: -1, B: 1})

R2 = Reaction('R2')
R2.name = 'Reaction 2'
R2.lower_bound = 0
R2.upper_bound = 1000
R2.add_metabolites({B: -1, C: 1})

R3 = Reaction('R3')
R3.name = 'Reaction 3'
R3.lower_bound = 0
R3.upper_bound = 1000
R3.add_metabolites({C: -1, D: 1, ATP: -1, ADP: 1})

R4 = Reaction('R4')
R4.name = 'Reaction 4'
R4.lower_bound = 0
R4.upper_bound = 1000
R4.add_metabolites({D: -1, E: 1})

R5 = Reaction('R5')
R5.name = 'ATP Synthase'
R5.lower_bound = 0
R5.upper_bound = 1000
R5.add_metabolites({E: -1, A: 1, ADP: -1, ATP: 1})

# Add reactions to the model
model.add_reactions([R1, R2, R3, R4, R5])

# Add a demand reaction for ATP
demand_ATP = Reaction('DM_atp_c')
demand_ATP.lower_bound = 1  # Force some demand for ATP to ensure flux through the cycle
demand_ATP.upper_bound = 1000
demand_ATP.add_metabolites({ATP: -1})

model.add_reactions([demand_ATP])

# Set objective to ATP production
model.objective = 'R5'
