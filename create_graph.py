import collections
import pandas as pd
from enum import Enum
import math
# RUNING FROM APP
from . import graphml_parser as gp
from . import txt_parser as tp
# RUNING DIRECTLY
#import graphml_parser as gp
#import txt_parser as tp
import importlib
import copy
import re


import os
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(BASE_DIR, "assets")





importlib.reload(gp)
importlib.reload(tp)


# returns a set of precursors and set of final products AND returns a set of precursors and final products
def identify_precs_and_final_prod(reactions): # TODO TODO TODO
    """
    Identifies, which of the molecules are the initial reactants
    (there is no reaction leading TO them)
    and which of them are the final products
    (there is no reacction leading FROM them)
    Args:
        ractions: List of reaction objects
    Returns:
        Set: Set molecules (initial reactants)
        Set: Set of molecules (final products)
        Set: Set of molecules (all molecules)
    """
    precursors = set()
    final_products = set()
    all_mols = set()
    # adding all the reactants into the set of precursors
    for reaction in reactions:
        precursors.add(reaction.source)
        final_products.add(reaction.target)
        
    all_mols = precursors.union(final_products)
       
    for reaction in reactions:
        precursors.discard(reaction.target)
        final_products.discard(reaction.source)
    
    for precursor in precursors:
        precursor.is_initial_reactant = True
    for final_prod in final_products:
        final_prod.is_final_product = True
    #print("PRINTING PRECURSORS:")
    #print(precursors)
    #print("PRINTING FINAL PRODUCTS:")
    #print(final_products)
    return precursors, final_products, all_mols
    
def calculate_rel_error(value, stderr):
    return math.inf if value == 0 else stderr/value

def calculate_stat_relevance(reactions):
    RELATIVE_ERROR_TRESHOLD = 100.0
    for react in reactions:
        rel_error = calculate_rel_error(react.value, react.stderr)
        if rel_error > RELATIVE_ERROR_TRESHOLD:
            react.stat_rel = gp.statistical_relevance.LOW
        else:
            react.stat_rel = gp.statistical_relevance.HIGH
            
def print_relevant(reactions):
    #print("PRINTING THE STATISTICALY RELEVANT DATA")
    for react in reactions:
        doPrint = True
        for product in react.products:
            if "Sink" in product:
                doPrint = False
        if doPrint:
            if react.stat_rel == gp.statistical_relevance.HIGH:
                pass
                
def print_all(reactions):
    print("PRINTING ALL THE REACTIONS")
    for react in reactions:
        print(react)
             
def prov_sort(reactions):
    sorted_reactions = sorted(reactions, key=lambda react: react.value, reverse=True)
    print("PRINTING THE SORTED LIST OF REACTIONS")
    print(sorted_reactions)
    return sorted_reactions   
                   
class Pruning:
    """
    Represents an instance of the pruning algorithm
    Attributes:
        reactions: List of reaction objects
        first_reactants: List of molecules, which are not (intermediate) products
        final_products: List of molecules, which are not (initial) reactants
        ...TODO
    """
    def __init__(self, 
                 reactions, 
                 first_reactants, 
                 final_products, 
                 molecules, 
                 remove_irrelevant_pre=False,
                 basic_pruning=False,
                 connecting=False,
                 adding_beyond_treshold=False,
                 threshold=0.75,
                 proportion=False):
        """
        Initializes the pruning algorithm class.
        Args:
            ractions: List of reaction objects
            first_reactants: List of molecules, which are not (intermediate) products
            final_products: List of molecules, which are not (initial) reactants
            ...TODO
        Returns:
            str: A greeting message.
        """
        self.reactions = copy.deepcopy(reactions)
        self.first_reactants = copy.deepcopy(first_reactants)
        self.final_products = copy.deepcopy(final_products)
        self.threshold = threshold
        self.proportion = proportion
        self.all_molecules = copy.deepcopy(molecules)
        self.basic_pruning = basic_pruning
        self.connecting = connecting
        self.adding_beyond_treshold = adding_beyond_treshold
        self.remove_irrelevant_pre = remove_irrelevant_pre
        self.removed_irr_pre = False
        self.basic_pruning_done = False
        self.connecting_done = False
        self.adding_beyond_treshold_done = False
        self.sorted = False
        self.pruned_reactions = set()
        self.mode = None
        self.sorted_reactions = list()
              
    def sort(self):
        """
        Sorts the attribute 'reactions' in the descending order according to the value
        of each reaction.
        """
        #sorted_reactions = sorted(self.reactions, key=lambda react: react.value, reverse=True)
        
        # SORTING NEWLY ACCORDING TO THE ABSOLUTE VALUE
        sorted_reactions = sorted(self.reactions, key=lambda react: abs(react.value), reverse=True)
        self.sorted_reactions = sorted_reactions
        self.sorted = True
        return sorted_reactions
        
        
    # maybe totally WRONG?
    def check_relevance(self, node, react1, react2):
        for react in self.sorted_reactions:
            if react == react1 or react == react2:
                continue
            elif react.source.name == node.name or react.target.name == node.name:
                # NASTAVENI THRESHOLDU NA JEDNU TISICINU
                if react.value > 0.001:
                    return True
        return False
    
    def remove_irrelevant_nodes_and_edges(self):
        
        
        print(f"Number of sorted reactions AFTER SORTING: {len(self.sorted_reactions)}")
        print(f"Number of pruned reactions AFTER SORTING: {len(self.pruned_reactions)}")
        #matches = []
        #self.sort()
        #print(f"Number of sorted reactions BEFORE removal: {len(self.sorted_reactions)}")
        #print(f"Number of pruned reactions BEFORE removal: {len(self.pruned_reactions)}")
        for react in self.sorted_reactions.copy():
            if react.source.name.endswith(".pre"):
                expected_target = react.source.name[:-4]
                if react.target.name == expected_target:
                    for react2 in self.sorted_reactions.copy():
                        if (#react2 != react \
                            react2.source.name == react.target.name \
                            and react2.target.name == "Sink" \
                            and abs(react.value - react2.value) < 0.5):
                            
                            # if the react.target is relevant, it is not to be removed from the mix
                            if self.check_relevance(react.target, react, react2):
                                print("----------------------------------------------------")
                                print("WE ARE NOT REMOVING -- MOLECULE IS TOO IMPORTANT!!!!")
                                print(f"Molecule in question: " + str(react.target.name))
                                print("----------------------------------------------------")
                                continue
                            
                            if (react.reaction_path == "polar" \
                                or react.reaction_path == "tracer_path" \
                                or react.reaction_path == "polar_path"\
                                or react2.reaction_path == "polar" \
                                or react2.reaction_path == "tracer_path" \
                                or react2.reaction_path == "polar_path"\
                                    ):
                                print("----------------------------------------------------")
                                print("WE ARE NOT REMOVING -- reaction IS POLAR!!!!")
                                #print(f"Molecule in question: " + str(react))
                                print("----------------------------------------------------")
                                continue
                            
                            if react in self.sorted_reactions:
                                print(f"REMOVING: " + react.equation + " val: " + str(react.value))
                                self.sorted_reactions.remove(react)
                            if react in self.pruned_reactions:
                                self.pruned_reactions.remove(react)
                            
                            if react2 in self.sorted_reactions:
                                print(f"REMOVING: " + react2.equation + " val: " + str(react2.value))
                                self.sorted_reactions.remove(react2)
                            if react2 in self.pruned_reactions:
                                self.pruned_reactions.remove(react2)
                            #print("Removing reaction:")
                            #print(react2)

                            #print("Removing nodes:")
                            #print(react.source)
                            #print(react.target)
                            self.first_reactants.remove(react.source)
                            self.all_molecules.remove(react.source)
                            self.all_molecules.remove(react.target)
                            for react3 in self.sorted_reactions.copy():
                                if react3.source.name == react.target.name or react3.target.name == react.target.name:
                                    if react3.source.name == react.target.name:
                                        print("the reac.target is source ")
                                    if react3.target.name == react.target.name:
                                        print("the reac.target is target ")
                                        
                                        
                                        
                                        
                                    if react3 in self.sorted_reactions:
                                        print(f"REMOVING related eq: " + react3.equation + " val: " + str(react3.value))
                                        self.sorted_reactions.remove(react3)
                                    if react3 in self.pruned_reactions:
                                        self.pruned_reactions.remove(react3)
        self.removed_irr_pre = True
        #print(f"Number of sorted reactions AFTER removal: {len(self.sorted_reactions)}")   
        #print(f"Number of pruned reactions AFTER removal: {len(self.pruned_reactions)}")
        
    def keep_polars_and_tracers(self):
        i = 0
        for react in self.reactions:
            if react.reaction_path == "polar" \
                or react.reaction_path == "tracer_path" \
                or react.reaction_path == "polar_path":
                i += 1
                print(f"POLAR REACTION KEPT: " + str(react.equation))
                
                self.pruned_reactions.add(react)
        print(f"Number of polar reactions: {i}")
        
                
    # tady dat Kruskala    
    def prune(self):
        """
        Conducts pruning - a distant 'variant' of Kruskal algorithm.
        Returns:
            set: set of reaction objects, which represent the 
            reactions, which were decided to be kept
        """
        if self.basic_pruning_done:
            print("Basic pruning has already been conducted!")
            return 
        #sorted_reactions = self.sort()
        #print("PRINTING THE SORTED REACTION >>>>>>>>>>>>>>>>>>>>>>>>")
        #for react in sorted_reactions:
        #    print(f">>>{react.equation, react.value, react.source, react.target} <<<\n")
        #    if (react.target.name == 'Pal.out'):
        #        print("TADY BYL!!!")
        #print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        
        # COMMENT IF ALREADY SORTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #self.sort()
        sorted_reactions = self.sorted_reactions

        
        pruned_reactions = set()
        
        last_reaction_value = 0
        for reaction in sorted_reactions:
            

            if (reaction.source.is_initial_reactant == False):  # reaktant není počáteční - musí vést do něj i z něj
                if (reaction.source.hasSource == False):
                    reaction.source.hasSource = True
                    pruned_reactions.add(reaction)
                    #try: 
                    #    sorted_reactions.remove(reaction)
                    #except ValueError as e:
                    #    pass
                        #print("Reaction has already been removed.")
                    #last_reaction_value = reaction.value
            if (reaction.source.hasTarget == False):  # ať už reaktant je či není počáteční, musí z něj vést hrana
                reaction.source.hasTarget = True
                pruned_reactions.add(reaction)
                #try: 
                #    sorted_reactions.remove(reaction)
                #except ValueError as e:
                #    pass
                    #print("Reaction has already been removed.")
                #last_reaction_value = reaction.value

            if (reaction.target.is_final_product == False):
                if (reaction.target.hasTarget == False):
                    reaction.target.hasTarget = True
                    pruned_reactions.add(reaction)
                    #try: 
                    #    sorted_reactions.remove(reaction)
                    #except ValueError as e:
                    #    pass
                        #print("Reaction has already been removed.")
                    #last_reaction_value = reaction.value   
            if (reaction.target.hasSource == False):
                reaction.target.hasSource = True
                pruned_reactions.add(reaction)
                #try: 
                #    sorted_reactions.remove(reaction)
                #except ValueError as e:
                #    pass
                    #print("Reaction has already been removed.")
                #last_reaction_value = reaction.value

        self.pruned_reactions = pruned_reactions
        self.mode = "Basic_pruning"
        return pruned_reactions#, sorted_reactions
    

    def doBFS(self, molecule):
        """
        Does a variation of Breadth-First Search (BFS) - all reached molecules
        have their attribute 'isPartOfConnectedGraph' set to 'True'
        Args:
            molecule: Molecule object, start of BFS
        Returns:
            set: set of visited molecules
        """
        molecules_visited = set()
        molecules_queue = []
        molecules_queue.append(molecule)
        #print(f"DOING BFS FOR MOLECULE: {molecule}")
        molecule.isPartOfConnectedGraph = True
        while (len(molecules_queue) != 0):
            mol = molecules_queue.pop()
            molecules_visited.add(mol)
            for reaction in self.pruned_reactions:
                if (reaction.source == mol and reaction.target.isPartOfConnectedGraph == False):
                    reaction.target.isPartOfConnectedGraph = True
                    molecules_queue.append(reaction.target)
                if (reaction.target == mol and reaction.source.isPartOfConnectedGraph == False):
                    reaction.source.isPartOfConnectedGraph = True
                    molecules_queue.append(reaction.source)
        return molecules_visited
    
    def connect(self):
        """
        Conducts the act of 'connecting' the components of the graph.
        Returns:
            molecule: the molecule, which is part of a yet undiscovered component is
            returned. The molecule is part of the reaction (source or target) which
            'connects' the two components. If no reaction is to be added, the 
            function return 'None'.
            
        """
        # TADY TED MENIM - misto:
        #sorted_reactions = self.sort()
        # BERU JIZ SERTONENE
        sorted_reactions = self.sorted_reactions
        for reaction in sorted_reactions:
            if reaction.source.isPartOfConnectedGraph == True \
                and reaction.target.isPartOfConnectedGraph == False:
                    self.pruned_reactions.add(reaction)
                    #print(f"RETURNING FROM CONNECT (target): {reaction.target}")
                    return reaction.target
            if reaction.source.isPartOfConnectedGraph == False \
                and reaction.target.isPartOfConnectedGraph == True:
                    self.pruned_reactions.add(reaction)
                    #print(f"RETURNING FROM CONNECT (source): {reaction.source}")
                    return reaction.source  
        return None
          
    def ensure_connectivity(self):
        """
        Ensures that the graph is a connected graph.
        - first molecule is taken, for which BFS is conducted
        - through BFS, all the molecules, which are 'reachable' from
        the first molecule are marked as seen (attribute isPartOfConnectedGraph)
        - then a check is conducted, whether all the molecules are 
        per of the connected graph
        - if some are not, one of them is chosen by the function connect,
        and through BFS, all the other molecules in the component are marked as 
        connected (attribut isPartOfConnectedGraph)
        """
        if self.connecting_done:
            print("Connectivity has already been ensured!")
        list_of_molecules = self.all_molecules.copy()
        #print(f"CONNECTIVITY BEING DONE NOW, list len: {len(list_of_molecules)}")
        set_of_all_molecules = set(list_of_molecules)
        molecule = list_of_molecules.pop() 
        molecules_visited = set()
        i = 0
        search = True
        self.mode = "Connected_pruning" 
        #print(f"CONNECTIVITY BEING FINISHED0, list len: {len(list_of_molecules)}")
        while (search):
            i += 1
            molecules_newly_visited = self.doBFS(molecule)
            molecules_visited.union(molecules_newly_visited)
            if (molecules_visited == set_of_all_molecules):
                #print("Conectivity has been finished")
                search = False
            molecule = self.connect()
            if molecule is None:
                #print("We are connected!!!")
                #print(f"CONNECTIVITY BEING FINISHED2, set len: {len(set_of_all_molecules)}")
                #print(f"CONNECTIVITY BEING FINISHED2, list len: {len(molecules_visited)}")

                break
            #molecule = set_of_all_molecules.difference(molecules_visited)
            if i > 1000:
                print("Zacyklení ve while!!!")
                break
            
    def addBeyondTreshold(self, threshold=0.75):
        """
        Adds more potentially relevant reactions (edges) to the graph.
        Args:
            threshold: decides the quantile, which determines which
            reactions (according to their reaction value) are to be 
            added. If self.proportion == True, it adds reactions until the 
             percentage of pruned_reactions matches the threshold
        """
        if self.proportion == True:
            for react in self.sorted_reactions:
                if len(self.pruned_reactions) < len(self.sorted_reactions) * self.threshold:
                    self.pruned_reactions.add(react)
                else:
                    return
        if self.adding_beyond_treshold_done:
            print("You have already added reactions beyond certain treshold. \n\
                If you want to add more reactions, please initiate \
                new instance of algorithm with a different threshold.")
        self.treshold = threshold
        values = []
        for react in self.sorted_reactions:
            values.append(react.value)
            
        data_frame = pd.DataFrame(values, columns=["value"])
        #quantile = float(data_frame.iloc[0].quantile(threshold))
        #print(data_frame.iloc[1])
        quantile = float(data_frame["value"].quantile(threshold))
        for react in self.sorted_reactions:
            if react.value > quantile:
                self.pruned_reactions.add(react)
        self.mode = "threshold_" + str(int(threshold * 100))
        self.adding_beyond_treshold_done = True
        
    def print_statement(self):
        num_kept = len(self.reactions)
        if self.pruned_reactions != None:
            num_kept = len(self.pruned_reactions)
        percentage = num_kept/len(self.reactions) * 100
        print(f"Number of reactions kept: {num_kept}, which is {percentage:.2f}% of all reactions.")
        
    def run(self):

        print(f"Number of sorted reactions START: {len(self.sorted_reactions)}")   
        print(f"Number of pruned reactions START: {len(self.pruned_reactions)}")   
        
        if self.sorted:
            print("Already sorted.")
        else:
            self.sort()
        
         # FUNKCE ODSTRANUJICI IRELEVANTNI NODES A EDGES
        # ----------------------------------------------------
        if self.remove_irrelevant_pre:
            if self.removed_irr_pre:
                print("You have already eliminated reactions with x.pre -1-> x -1-> Sink")
            else:
                self.remove_irrelevant_nodes_and_edges()
        # ----------------------------------------------------
        print(f"Number of sorted reactions AFTER REMOVING IRRELEVANT: {len(self.sorted_reactions)}")   
        print(f"Number of pruned reactions AFTER REMOVING IRRELEVANT: {len(self.pruned_reactions)}")   
        
        
        self.keep_polars_and_tracers()
        
        print(f"Number of sorted reactions AFTER KEEPING POLARS AND TRACERS: {len(self.sorted_reactions)}")   
        print(f"Number of pruned reactions AFTER KEEPING POLARS AND TRACERS: {len(self.pruned_reactions)}")   
        
        
        if self.basic_pruning:
            if self.basic_pruning_done:
                print("Basic pruning has already been conducted!")
            else:
                self.prune()
                self.basic_pruning_done = True
                print("Basic pruning has been conducted.")
                self.print_statement()
                
       
        print(f"Number of sorted reactions AFTER BASIC PRUNING: {len(self.sorted_reactions)}")   
        print(f"Number of pruned reactions AFTER BASIC PRUNING: {len(self.pruned_reactions)}")   
        
        
        
        if self.connecting:
            if self.connecting_done:
                print("Connectivity has already been ensured!")
            else:
                self.ensure_connectivity()
                self.connecting_done = True
                print("Connectivity has been ensured.")
                self.print_statement()
                
        print(f"Number of sorted reactions AFTER CONNECTING: {len(self.sorted_reactions)}")   
        print(f"Number of pruned reactions AFTER CONNECTING: {len(self.pruned_reactions)}")           
                
                
        if self.adding_beyond_treshold:
            if self.adding_beyond_treshold_done:
                print("You have already added reactions beyond certain treshold. \n\
                    If you want to add more reactions, please initiate \
                    new instance of algorithm with a different threshold.")
            else:
                self.addBeyondTreshold(self.threshold)
                self.adding_beyond_treshold_done = True
                print(f"All reactions beyond threshold {self.threshold} have been added.")
                self.print_statement()
        
        print(f"Number of sorted reactions AFTER ADDING BEYOND THRESHOLD: {len(self.sorted_reactions)}")   
        print(f"Number of pruned reactions AFTER ADDING BEYOND THRESHOLD: {len(self.pruned_reactions)}")   
        
        self.keep_polars_and_tracers()  
                

# VISUALISATION CLASS

class Visualisation:
    """
    Represents an instance of the pruning algorithm
    Attributes:
        reactions: List of reaction objects
        first_reactants: List of molecules, which are not (intermediate) products
        final_products: List of molecules, which are not (initial) reactants
        ...TODO
    """
    def __init__(self, 
                 reactions, 
                 first_reactants, 
                 final_products, 
                 molecules, 
                 basic_pruning=False,
                 connecting=False,
                 adding_beyond_treshold=False,
                 init_num_reactions=10,
                 proportion=False):
        """
        Initializes the pruning algorithm class.
        Args:
            ractions: List of reaction objects
            first_reactants: List of molecules, which are not (intermediate) products
            final_products: List of molecules, which are not (initial) reactants
            ...TODO
        Returns:
            str: A greeting message.
        """
        self.reactions = copy.deepcopy(reactions)
        self.first_reactants = copy.deepcopy(first_reactants)
        self.final_products = copy.deepcopy(final_products)
        self.init_num_reactions = init_num_reactions
        self.proportion = proportion
        self.all_molecules = copy.deepcopy(molecules)
        self.basic_pruning = basic_pruning
        self.connecting = connecting
        self.adding_beyond_treshold = adding_beyond_treshold
        self.basic_pruning_done = False
        self.connecting_done = False
        self.adding_beyond_treshold_done = False
        self.pruned_reactions = set()
        self.mode = "DEBUGGING"
        self.sorted_reactions = list()
    
    
    def sort(self):
        """
        Sorts the attribute 'reactions' in the descending order according to the value
        of each reaction.
        """
        sorted_reactions = sorted(self.reactions, key=lambda react: react.value, reverse=True)
        self.sorted_reactions = sorted_reactions
        return sorted_reactions
        
    # MUST BE DONE BEFORE BFS
    def add_edges(self):
        """
        Adding first x (init_num_reactions) reactions with the highest value (react.value)
        """
        #sorted_reactions = self.sort()
        sorted_reactions = self.sorted_reactions
        #print("PRINTING ADDED REACTIONS")
        for i in range(self.init_num_reactions):
            self.pruned_reactions.add(sorted_reactions[i])
            sorted_reactions[i].source.hasTarget = True
            sorted_reactions[i].target.hasSource = True
            #print(sorted_reactions[i])
            
    # MUST BE DONE BEFORE BFS
    def include_specific_node(self, node_name):
        """
        Adding edges, which ought to be included.
        """
        for react in self.sorted_reactions: # projizdim reakcemi
            
            if react.source.name == node_name:  # najdu reakci, kde je zdrojem hledana molekula
                # NEMUZE SE STAT, ze by dana molekula nemela vychazejici reakci (protoze je source v react)
                # zaroven je to nejvetsi takova reakce (protoze ji potkame jako prvni)
                if (react.source.hasTarget == False): # nema cil? tak pridam tuto reakci
                    self.pruned_reactions.add(react)
                    react.source.hasTarget = True
                    
            if react.target.name == node_name:
                # NEMUZE SE STAT, ze by dana molekula nemela prichazejici reakci (protoze je target v react)
                # zaroven je to nejvetsi takova reakce (protoze ji potkame jako prvni)
                if (react.target.hasSource == False):
                    self.pruned_reactions.add(react)
                    react.target.hasSource = True
        print(f"The following node should be included: {node_name}")
    
    def include_specific_nodes(self, nodes):
        for single_node in nodes:
            self.include_specific_node(single_node)
        
            
    def add_bfs(self):
        reactions_added = set()
        reactions_queue = []
        reactions_queue.extend(self.pruned_reactions)
        while(len(reactions_queue) != 0):
            react = reactions_queue.pop()
            reactions_added.add(react)
            if (react.source.is_initial_reactant == False):
                if react.source.hasSource == False:
                    for pot_react in self.sorted_reactions:
                        if pot_react.target == react.source:
                            react.source.hasSource = True
                            if pot_react not in reactions_added:
                                reactions_queue.append(pot_react)
                            break
            if (react.target.is_final_product == False):
                if react.target.hasTarget == False:
                    for pot_react in self.sorted_reactions:
                        if pot_react.source == react.target:
                            react.target.hasTarget = True
                            if pot_react not in reactions_added:
                                reactions_queue.append(pot_react)
                            break
        self.pruned_reactions = reactions_added
        return reactions_added
    
    # DLE JMENA!!
    def unite_molecules(self, molecules_to_unite, new_name):
        new_mol = gp.molecule(new_name)
        for molecule in molecules_to_unite:
            for react in self.pruned_reactions:
                if react.source.name == molecule: 
                    #print("RENAMING")
                    if new_mol.id == 0:
                        new_mol.id = react.source.id
                    #print("PRINTING THE RENAMED MOLECULE and the reaction:")
                    new_mol.addOldNameId(react.source.name, react.source.id)
                    react.source = new_mol 
                    pattern = r'\b' + re.escape(molecule) + r'\b'
                    react.equation = re.sub(pattern, new_name, react.equation)
                    
                    #react.equation = react.equation.replace(molecule, new_name)   
                    #print(react.source)
                    #print(react)
                if react.target.name == molecule: 
                    #print("RENAMING")
                    if new_mol.id == 0:
                        new_mol.id = react.target.id
                    #print("PRINTING THE RENAMED MOLECULE and the reaction:")
                    new_mol.addOldNameId(react.target.name, react.target.id)
                    react.target = new_mol
                    
                    pattern = r'\b' + re.escape(molecule) + r'\b'
                    react.equation = re.sub(pattern, new_name, react.equation)
                    #react.equation = react.equation.replace(molecule, new_name)   
                    #print(react.target)
                    #print(react)
            
    def remove_self_loops(self):
        for react in self.pruned_reactions.copy():
            if react.source == react.target:
                self.pruned_reactions.remove(react)
                
    def unite_reactions(self):
        for react1 in self.pruned_reactions.copy():
            for react2 in self.pruned_reactions.copy():
                if react1 != react2:
                    if react1.source == react2.source and react1.target == react2.target:
                        react1.value = react1.value + react2.value
                        self.pruned_reactions.remove(react2)
                    if react1.source == react2.target and react1.target == react2.source:
                        if react1.value > react2.value:
                            react1.value = react1.value - react2.value
                            self.pruned_reactions.remove(react2)
                        else:
                            react2.value = react2.value - react1.value
                            self.pruned_reactions.remove(react1)
                            
        
    
                
                    
                
    
    
        
    
        
        
# RUNING THE RESPECTIVE PARTS OF THE ALGORITHM:
# the following code is connected to the "pruning" algorithm (class prune)     


def prepare_mols_reacs(data_address_input):
    INPUT_PATH_GRAPHML = os.path.join(ASSETS_DIR, "input", data_address_input)
    molecules, reactions, node_id_to_mol = gp.data_parser(INPUT_PATH_GRAPHML)
    #molecules, reactions, node_id_to_mol = gp.data_parser(data_address_input)
    first_reactants, final_products, all_mols = identify_precs_and_final_prod(reactions)
    return all_mols, reactions, first_reactants, final_products, node_id_to_mol

   
def save_result_graphml(INPUT_PATH_GRAPHML, OUTPUT_PATH, first_algorithm, remove_nodes=False, name="output"):
    INPUT_PATH_GRAPHML = os.path.join(ASSETS_DIR, "input", INPUT_PATH_GRAPHML)
    #molecules, reactions, node_id_to_mol = gp.data_parser(INPUT_PATH_GRAPHML)
    
    
    #OUTPUT_PATH_GRAPHML = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".graphml"
    # TOTO naposledy zakomentovano
    #OUTPUT_PATH_GRAPHML = f"" + OUTPUT_PATH + "/" + name + "_T" + str(int(first_algorithm.threshold * 100)) + ".graphml"
    
    OUTPUT_PATH = os.path.join(ASSETS_DIR, "output", OUTPUT_PATH + "_T" + str(int(first_algorithm.threshold * 100)) + ".graphml")
    gp.to_graphml(INPUT_PATH_GRAPHML, OUTPUT_PATH, first_algorithm.pruned_reactions, remove_nodes)
    #gp.to_graphml(INPUT_PATH_GRAPHML, OUTPUT_PATH_GRAPHML, first_algorithm.pruned_reactions, remove_nodes)
    
def save_result_graphml_vis(INPUT_PATH_GRAPHML, OUTPUT_PATH, first_algorithm, name="output"):
    INPUT_PATH_GRAPHML = os.path.join(ASSETS_DIR, "input", INPUT_PATH_GRAPHML)
    OUTPUT_PATH_GRAPHML = os.path.join(ASSETS_DIR, "output", OUTPUT_PATH)
    
    #OUTPUT_PATH_GRAPHML = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".graphml"
    #OUTPUT_PATH_GRAPHML = f"" + OUTPUT_PATH + "/" + name + ".graphml"
    gp.to_graphml_visualize(INPUT_PATH_GRAPHML, OUTPUT_PATH_GRAPHML, first_algorithm.pruned_reactions)
    
def save_result_txt(INPUT_PATH_TXT, OUTPUT_PATH, first_algorithm, name="output"):
    #OUTPUT_PATH_TXT = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".txt"
    OUTPUT_PATH_TXT = os.path.join(ASSETS_DIR, "output", OUTPUT_PATH + "_T" + str(int(first_algorithm.threshold * 100)) + ".txt")
    INPUT_PATH_TXT = os.path.join(ASSETS_DIR, "input", INPUT_PATH_TXT)
    #OUTPUT_PATH_TXT = f"" + OUTPUT_PATH + "/" + name + "_T" + str(int(first_algorithm.threshold * 100)) + ".txt"
    tp.run_the_script(first_algorithm.pruned_reactions, first_algorithm.reactions, INPUT_PATH_TXT, OUTPUT_PATH_TXT)
    

    
def run_script(INPUT_PATH_GRAPHML, INPUT_PATH_TXT, OUTPUT_PATH, name="output"):
    # PREPARING THE ALGORITHM
    molecules, reactions, node_id_to_mol = gp.data_parser(INPUT_PATH_GRAPHML)
    calculate_stat_relevance(reactions) # can be omitted
    first_reactants, final_products, all_mols = identify_precs_and_final_prod(reactions)
    first_algorithm = Pruning(reactions, 
                              first_reactants, 
                              final_products, 
                              all_mols,
                              basic_pruning=True,
                              connecting=True,
                              adding_beyond_treshold=True,
                              threshold=0.75
                              )
    first_algorithm.run()
    
    # CONDUCTING THE PRUNING
    #first_algorithm.prune()
    #pruned_reactions = first_algorithm.pruned_reactions
    #print(f"Number of pruned reactions after PRUNING: {len(pruned_reactions)}")
    #
    ## CONNECTING 
    #first_algorithm.ensure_connectivity()
    #pruned_reactions = first_algorithm.pruned_reactions
    #print(f"Number of pruned reactions after PRUNING + CONNECTING: {len(pruned_reactions)}")
    #
    ## ADDING BEYOND THRESHOLD
    #first_algorithm.addBeyondTreshold(threshold=0.50)   # deciding the threshold
    #pruned_reactions = first_algorithm.pruned_reactions
    #print(f"Number of pruned reactions after PRUNING + CONNECTING + ADDING QUANTILES: {len(pruned_reactions)}")
    
    ## SAVING TO .GRAPHML FILE
    #save_result_graphml(OUTPUT_PATH, first_algorithm, name="testing", keep_all_mols=True, node_id_to_mol=node_id_to_mol)
    #OUTPUT_PATH_GRAPHML = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".graphml"
    #gp.to_graphml(INPUT_PATH_GRAPHML, OUTPUT_PATH_GRAPHML, pruned_reactions)
    #
    ## SAVING TO .TXT FILE
    #save_result_txt(INPUT_PATH_TXT, OUTPUT_PATH, first_algorithm, name="testing")
    
    
    
    
    #OUTPUT_PATH_TXT = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".txt"
    #ax.run_the_script(pruned_reactions, first_algorithm.reactions, INPUT_PATH_TXT, OUTPUT_PATH_TXT)

# SECOND ALGORITHM

def save_res(INPUT_PATH_GRAPHML, OUTPUT_PATH, second_algorithm, name="NONAME"):
    save_result_graphml_vis(INPUT_PATH_GRAPHML, OUTPUT_PATH, second_algorithm, name=name)
    

def run_sec_algo():
    # path to the source graphml file - the graph to be pruned
    INPUT_PATH_GRAPHML = "./metabolomic_optimization/assets/input/PalPaoSteOle_regular.graphml"
    # path to the original reactions - reactions to be selected
    INPUT_PATH_TXT = "./metabolomic_optimization/assets/input/PalPaoSteOle_regular_new.txt"
    # path to the output, where both the resulting .graphml file and .txt file is stored
    OUTPUT_PATH = "./metabolomic_optimization/assets/output"


    molecules, reactions, node_id_to_mol = gp.data_parser(INPUT_PATH_GRAPHML)
    calculate_stat_relevance(reactions) # can be omitted
    first_reactants, final_products, all_mols = identify_precs_and_final_prod(reactions)

    second_algorithm = Visualisation(
        reactions, 
        first_reactants, 
        final_products, 
        all_mols,
        basic_pruning=True,
        connecting=True,
        adding_beyond_treshold=True,
        init_num_reactions=10
    )

    second_algorithm.sort()
    second_algorithm.add_edges()
    added_reactions = second_algorithm.add_bfs()
    #list_of_mols_unite = ["TG123StePaoPal.pre", "TG123StePaoPal"]

    list_of_mols_unite = [
        "GAP", 
        "Pyr.c", 
        "Fum.m", 
        "Mal.m", 
        "Pyr.m", 
        "Ac.m", 
        "Oac.m", 
        "Cit.m", 
        "Acon", 
        "Cit.c"
    ]

    list_of_mols_unite2 = [
        "Glc.lab",
        "Glc",
        "DHAP",
        "G3P.New",
        "G3P"
    ]

    list_of_mols_unite3 = [
        "Pao.P2",
        "TG123PaoStePao"
    ]

    list_of_mols_unite4 = [
        "Ac",
        "Pal.Pre"
    ]

    # TADY KOMENTOVAT/NEKOMENTOVAT PRO SDRUZOVANI MOLEKUL
    #second_algorithm.unite_molecules(list_of_mols_unite, "CITRATOVY CYKLUS")
    second_algorithm.unite_molecules(list_of_mols_unite2, "GLYKOLYSA?")
    ##second_algorithm.unite_molecules(list_of_mols_unite3, "!!!TEST!!!")
    #second_algorithm.unite_molecules(list_of_mols_unite4, "OOOOOOOOOOOOOOOOOOOOOOOOO")
    second_algorithm.remove_self_loops()
    #second_algorithm.unite_reactions()

    print_all(second_algorithm.pruned_reactions)

    #save_result_txt(INPUT_PATH_TXT, OUTPUT_PATH, second_algorithm, name="UNITED")


    save_res(INPUT_PATH_GRAPHML, OUTPUT_PATH, second_algorithm, name="UNITED2")


# TADY RUZNE KOMENTOVAT a ZAKOMENTOVAT
#run_sec_algo()


# COMMENT THE FOLOOWING IN ORDER FOR THE SCRIPT TO WORK!
"""
# path to the source graphml file - the graph to be pruned
INPUT_PATH_GRAPHML = "./metabolomic_optimization/assets/input/threepath.graphml"
# path to the original reactions - reactions to be selected
INPUT_PATH_TXT = "./metabolomic_optimization/assets/input/PalPaoSteOle_regular.txt"
# path to the output, where both the resulting .graphml file and .txt file is stored
OUTPUT_PATH = "./metabolomic_optimization/assets/output"

run_script(INPUT_PATH_GRAPHML,INPUT_PATH_TXT, OUTPUT_PATH)

"""