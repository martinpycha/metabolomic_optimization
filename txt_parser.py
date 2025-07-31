import pandas as pd
import os
import re

def add_header(file_location):
    # Load the existing Excel file (no headers assumed)
    df = pd.read_excel(file_location, header=None)
    # Define the new header row
    new_header = ["ID", "Equation", "Value", "StdErr", "LB", "UB", "Pathway"]
    # Insert the new header row at the top
    df.columns = new_header  # assign as column names
    df.dropna(subset=['Equation'], inplace=True)
    # Save the new file (overwrite or save separately)
    df.to_excel("PalPaoSteOle_regular_1279_flux_v2_withHeader.xlsx", index=False)

#add_header("./metabolomic_optimization/PalPaoSteOle_regular_1279_flux_v2.xlsx")

def read_txt_reaction_file(file_location):
    reactions_dict = {}
    with open(file_location) as file:
        reactions = [line.rstrip() for line in file]        
    for original_react in reactions:
        react = re.sub(r'\s*\([^)]*\)', '', original_react)
        react = react.split(';')[0]
        react = react.strip("'\"").strip()
        reactions_dict[react] = original_react
    return reactions_dict

def keep_only_relevant(reactions_dict, pruned_reactions, all_reactions, output_file_path):
    reactions_to_write = list()
    already_added = set()
    #print("Writing the reactions:")
    for reaction, original_reaction in reactions_dict.items():
        #print(f"Reaction: {original_reaction}")
        for reaction_p in pruned_reactions:
            if reaction.strip() == reaction_p.equation.strip() and original_reaction not in already_added:
                reactions_to_write.append(original_reaction)
                already_added.add(original_reaction)
                
    #print("End of writing the reactions")
    #pruned_reactions_equations = set()
    #all_reactions_equations = set()
    #for r in pruned_reactions:
    #    if r.equation in pruned_reactions_equations:
    #        continue
    #    else:
    #        pruned_reactions_equations.add(r.equation)
    #    pruned_reactions_equations.add(r.equation)
    #for r in all_reactions:
    #    all_reactions_equations.add(r.equation)
    #print(f"{len(pruned_reactions_equations)}")
    reactions_dict_keys = set(reactions_dict.keys())
    with open(output_file_path, 'w') as file:
        for equation in reactions_to_write:
            file.write('%s\n' %equation)    
    print("Succesfully written the result into the file!")
    
def run_the_script(pruned_reactions, all_reactions, input_file_path, output_file_path):
    #reactions_dict = read_txt_reaction_file("./metabolomic_optimization/PalPaoSteOle_regular.txt")
    reactions_dict = read_txt_reaction_file(input_file_path)
    keep_only_relevant(reactions_dict, pruned_reactions, all_reactions, output_file_path)