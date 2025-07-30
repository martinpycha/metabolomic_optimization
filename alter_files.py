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
        num_lines = len(reactions)
        
    set_of_reactions = set()
    for original_react in reactions:
        react = re.sub(r'\s*\([^)]*\)', '', original_react)
        react = react.split(';')[0]
        react = react.strip("'\"").strip()
        reactions_dict[react] = original_react
        set_of_reactions.add(react)
    return reactions_dict

def keep_only_relevant(reactions_dict, pruned_reactions, all_reactions, output_file_path):
    reactions_to_write = set()
    for reaction, original_reaction in reactions_dict.items():
        for reaction_p in pruned_reactions:
            if reaction.strip() == reaction_p.equation.strip():
                reactions_to_write.add(original_reaction)
    pruned_reactions_equations = set()
    all_reactions_equations = set()
    for r in pruned_reactions:
        if r.equation in pruned_reactions_equations:
            continue
        #    print("Already added.")
        #    print(f"   NOT adding the following reaction: {r.equation} with ID: {r.id} ")
        else:
            pruned_reactions_equations.add(r.equation)
        #    print(f"\nAdding the following reaction: {r.equation} with ID: {r.id}")
        pruned_reactions_equations.add(r.equation)
    for r in all_reactions:
        all_reactions_equations.add(r.equation)
    #print(f"\nLength of reactions in txt: {len(reactions_dict)}, \
    #        \nLength of all reactions from GRAPHML: {len(all_reactions)}, \
    #        \nLength of reactions pruned SET (ACCORDING TO ID): {len(pruned_reactions)}, \
    #        \nNumber of reactions pruned SET (ACCORDING TO EQUATION): {len(pruned_reactions_equations)}, \
    #        \nLength of reactions both in txt and in reactions pruned: {len(reactions_to_write)}. \n")
    pruned_reactions_equations = set()
    all_reactions_equations = set()
    for r in pruned_reactions:
        if r.equation in pruned_reactions_equations:
            continue
        #    print("Already added.")
        #    print(f"   NOT adding the following reaction: {r.equation} with ID: {r.id} ")
        else:
            pruned_reactions_equations.add(r.equation)
        #    print(f"\nAdding the following reaction: {r.equation} with ID: {r.id}")
        pruned_reactions_equations.add(r.equation)
    for r in all_reactions:
        all_reactions_equations.add(r.equation)
    print(f"{len(pruned_reactions_equations)}")
    reactions_dict_keys = set(reactions_dict.keys())
    #print(f"Number of reactions, which are in original txt but are not in the pruned reactions: {len(not_in_set)} \n \
    #        \n Number of reactions, which are in the pruned reaction, but are not in the reactions dictionary: {len(not_in_dict)} \n \
    #        \n Number of reactions both in Inca and original txt: {len(in_inca_and_txt)}")
    #print(f"COUNTING !!!! \n")
    #print(f"{len(pruned_reactions_equations)}")
    with open(output_file_path, 'w') as file:
        for equation in reactions_to_write:
            file.write('%s\n' %equation)    
    print("Succesfully written the result into the file!")
    
def run_the_script(pruned_reactions, all_reactions, input_file_path, output_file_path):
    #reactions_dict = read_txt_reaction_file("./metabolomic_optimization/PalPaoSteOle_regular.txt")
    reactions_dict = read_txt_reaction_file(input_file_path)
    keep_only_relevant(reactions_dict, pruned_reactions, all_reactions, output_file_path)