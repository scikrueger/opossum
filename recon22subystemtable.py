#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 20:22:52 2024

@author: flo
"""

import cobra
import pandas as pd


mod = cobra.io.read_sbml_model("/home/flo/Documents/resources/Recon22.xml")
mod = cobra.io.read_sbml_model("/home/flo/Documents/resources/Recon3DModel_301.xml")
mod.solver = "cplex"


# Create a list to store reaction data
data = []

for reaction in mod.reactions:
    # Extract 'NOTES' from the notes dictionary, handle the absence of 'NOTES'
    notes_content = reaction.notes.get('NOTES', 'No notes available')  # Default to 'No notes available' if 'NOTES' is not found

    # Collecting data for each reaction
    reaction_data = {
        "reaction": reaction.id,
        "name": reaction.name,
        "reaction eq": str(reaction.reaction),
        "subsystem": reaction.subsystem,
        "notes": notes_content  # Only extracting the 'NOTES' subcategory
    }
    data.append(reaction_data)
    
df = pd.DataFrame(data)

df.to_csv("~/Documents/resources/Recon22_metadata.csv", index = False)
df.to_csv("~/Documents/resources/Recon3D_metadata.csv", index = False)
