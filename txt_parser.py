import pandas as pd
import os
import re

def add_header(file_location):
    # Load the existing Excel file (no headers assumed)
    df = pd.read_excel(file_location, header=None)
    # Define the new header row
    #new_header = ["ID", "Equation", "Value", "StdErr", "LB", "UB", "Pathway"]
    new_header = ["Equation", "ID", "Value", "StdErr", "Pathway"]
    # Insert the new header row at the top
    df.columns = new_header  # assign as column names
    df.dropna(subset=['Equation'], inplace=True)
    # Save the new file (overwrite or save separately)
    
    # THERE WAS A MISTAKE - nets are to stay, exch are to be removed!
    #df = df[~df['ID'].astype(str).str.contains('net', case=False, na=False)]
    df = df[~df['ID'].astype(str).str.contains('exch', case=False, na=False)]
    if 'LB' not in df.columns:
        df['LB'] = ""
    if 'UB' not in df.columns:
        df['UB'] = ""
        
    desired_order = ["ID", "Equation", "Value", "StdErr", "LB", "UB", "Pathway"]
    df = df[desired_order]

    
    #df.to_excel("A4GALT_flux_withHeader.xlsx", index=False)
    df.to_excel("metabolomic_optimization/assets/input/PPSO_v2_header.xlsx", index=False)
    
#add_header("metabolomic_optimization/assets/input/PPSO_BEST.xlsx")
#add_header("metabolomic_optimization/PalPaoSteOle_regular_1279_flux_v2.xlsx")
    

    
def compare_two_txt_files(file1_location, file2_location):
    list1, list2 = list(), list()
    file1 = open(file1_location, 'r')
    for line in file1:
        list1.append(line)
    file2 = open(file2_location, 'r')
    for line in file2:
        list2.append(line)
    diff1 = [item for item in list1 if item not in list2]
    diff2 = [item for item in list2 if item not in list1]
    print("-------------------------------------------------------------------------")
    print(f"First file has {len(list1)} rows, Second file has {len(list2)} rows.")
    print("-------------------------------------------------------------------------")
    print(f"The following reactions are in the first file and not in the second one:")
    if len(diff1) == 0:
        print("All reactions in second file are also in the first line.")
    for react in diff1:
        print(f"React: {react}")
    #print(diff1)
    print("-------------------------------------------------------------------------")
    print(f"The following reactions are in the second file and not in the first one:")
    if len(diff2) == 0:
        print("All reactions in first file are also in the second line.")
    for react in diff2:
        print(f"React: {react}")
    #print(diff2)
    print("-------------------------------------------------------------------------")

#compare_two_txt_files("./metabolomic_optimization/assets/input/model_202r_output.txt", "./metabolomic_optimization/assets/output/A4GALT_threshold_0.txt")

# Revelation of the problem at the 40th percentil
#compare_two_txt_files("./metabolomic_optimization/assets/output/PalPaoSteOle_regular_new_threshold_40.txt", "./metabolomic_optimization/assets/input/PalPaoSteOle_regular_PROBLEM_threshold_40.txt")

import re
import pandas as pd

def clean_reaction(line):
    """
    Remove bracketed sections and extra spaces from a reaction line.
    Example:
    'MG1Pal (...) + Pal (...) -> DG12PalPal (...)'
    becomes
    'MG1Pal + Pal -> DG12PalPal'
    """
    # Remove everything inside parentheses
    cleaned = re.sub(r'\([^)]*\)', '', line)
    # Remove multiple spaces
    cleaned = re.sub(r'\s+', ' ', cleaned).strip()
    return cleaned

def compare_files(txt_file, excel_file, excel_column=1, output_file="comparison_results.txt"):
    # Read TXT file (long reactions with brackets)
    with open(txt_file, "r") as f:
        original_lines = [clean_reaction(line.strip()) for line in f if line.strip()]
    
    # Read Excel file (simplified reactions, assumed in first column by default)
    df = pd.read_excel(excel_file, header=None)
    simplified_lines = [str(x).strip() for x in df[excel_column].dropna().tolist()]
    
    # Compare
    extra_in_original = [r for r in original_lines if r not in simplified_lines]
    extra_in_simplified = [r for r in simplified_lines if r not in original_lines]
    matches = [r for r in simplified_lines if r in original_lines]
    
    # Write results to file
    with open(output_file, "w") as out:
        out.write(f"✅ Matches found: {len(matches)}\n\n")
        
        out.write("⚠️ Extra reactions in TXT (long) file:\n")
        for r in extra_in_original:
            out.write("  " + r + "\n")
        
        out.write("\n⚠️ Extra reactions in Excel (simplified) file:\n")
        for r in extra_in_simplified:
            out.write("  " + r + "\n")
    
    print(f"Comparison complete ✅ Results saved in {output_file}")

# Example usage:
#compare_files("./metabolomic_optimization/assets/input/PalPaoSteOle_regular_new.txt", "./metabolomic_optimization/assets/input/PPSO_regular_v2.xlsx", excel_column=0, output_file="results.txt")

#add_header("./metabolomic_optimization/A4GALT_flux.xlsx")

def read_txt_reaction_file(file_location):
    reactions_dict = {}
    with open(file_location) as file:
        reactions = [line.rstrip() for line in file]   
    print(f"Number of reactions parsed from .txt file: {len(reactions)}")     
    for original_react in reactions:
        react = re.sub(r'\s*\([^)]*\)', '', original_react)
        react = react.split(';')[0]
        react = react.strip("'\"").strip()
        
        if react in reactions_dict:
            print(f"FOR the reaction {react}, there are two reactions:")
            print(f"First one: {reactions_dict[react]}")
            print(f"Second one: {original_react}")
        reactions_dict[react] = original_react
    print(f"Number of reactions in reactions dict: {len(reactions_dict)}")

    return reactions_dict

def keep_only_relevant(reactions_dict, pruned_reactions, all_reactions, output_file_path):
    print(f"Length of pruned reactions: {len(pruned_reactions)}")
    
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
    print(f"Length of equation to write: {len(reactions_to_write)}")
    with open(output_file_path, 'w') as file:
        for equation in reactions_to_write:
            file.write('%s\n' %equation)    
    with open(output_file_path, 'r') as file:
        line_count = sum(1 for line in file)

    print(f'Total lines: {line_count}')
    print("Succesfully written the result into the file!")
    
def run_the_script(pruned_reactions, all_reactions, input_file_path, output_file_path):
    #reactions_dict = read_txt_reaction_file("./metabolomic_optimization/PalPaoSteOle_regular.txt")
    reactions_dict = read_txt_reaction_file(input_file_path)
    keep_only_relevant(reactions_dict, pruned_reactions, all_reactions, output_file_path)