import collections
import pandas as pd
from enum import Enum
import math
import graphml_parser as gp
import txt_parser as tp
import importlib
import copy

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
        self.basic_pruning_done = False
        self.connecting_done = False
        self.adding_beyond_treshold_done = False
        self.pruned_reactions = None
        self.mode = None
              
    def sort(self):
        """
        Sorts the attribute 'reactions' in the descending order according to the value
        of each reaction.
        """
        sorted_reactions = sorted(self.reactions, key=lambda react: react.value, reverse=True)
        return sorted_reactions
        
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
        sorted_reactions = self.sort()
        #print("PRINTING THE SORTED REACTION >>>>>>>>>>>>>>>>>>>>>>>>")
        #for react in sorted_reactions:
        #    print(f">>>{react.equation, react.value, react.source, react.target} <<<\n")
        #    if (react.target.name == 'Pal.out'):
        #        print("TADY BYL!!!")
        #print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        
        
        #sorted_reactions = self.reactions
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
        sorted_reactions = self.sort()
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
        set_of_all_molecules = set(list_of_molecules)
        molecule = list_of_molecules.pop() 
        molecules_visited = set()
        i = 0
        search = True
        self.mode = "Connected_pruning" 
        while (search):
            i += 1
            molecules_newly_visited = self.doBFS(molecule)
            molecules_visited.union(molecules_newly_visited)
            if (molecules_visited == set_of_all_molecules):
                search = False
            molecule = self.connect()
            if molecule is None:
                #print("We are connected!!!")
                break
            #molecule = set_of_all_molecules.difference(molecules_visited)
            if i > 1000:
                #print("Zacyklení ve while!!!")
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
            for react in self.reactions:
                if len(self.pruned_reactions) < len(self.reactions) * self.threshold:
                    self.pruned_reactions.add(react)
                else:
                    return
        if self.adding_beyond_treshold_done:
            print("You have already added reactions beyond certain treshold. \n\
                If you want to add more reactions, please initiate \
                new instance of algorithm with a different threshold.")
        self.treshold = threshold
        values = []
        for react in self.reactions:
            values.append(react.value)
            
        data_frame = pd.DataFrame(values, columns=["value"])
        #quantile = float(data_frame.iloc[0].quantile(threshold))
        #print(data_frame.iloc[1])
        quantile = float(data_frame["value"].quantile(threshold))
        for react in self.reactions:
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
        if self.basic_pruning:
            if self.basic_pruning_done:
                print("Basic pruning has already been conducted!")
            else:
                self.prune()
                self.basic_pruning_done = True
                print("Basic pruning has been conducted.")
                self.print_statement()
        if self.connecting:
            if self.connecting_done:
                print("Connectivity has already been ensured!")
            else:
                self.ensure_connectivity()
                self.connecting_done = True
                print("Connectivity has been ensured.")
                self.print_statement()
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
                
def prepare_mols_reacs(data_address_input):
    molecules, reactions, _ = gp.data_parser(data_address_input)
    first_reactants, final_products, all_mols = identify_precs_and_final_prod(reactions)
    return all_mols, reactions, first_reactants, final_products
   
def save_result_graphml(INPUT_PATH_GRAPHML, OUTPUT_PATH, first_algorithm, name="output"):
    OUTPUT_PATH_GRAPHML = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".graphml"
    gp.to_graphml(INPUT_PATH_GRAPHML, OUTPUT_PATH_GRAPHML, first_algorithm.pruned_reactions)
    
def save_result_txt(INPUT_PATH_TXT, OUTPUT_PATH, first_algorithm, name="output"):
    OUTPUT_PATH_TXT = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".txt"
    tp.run_the_script(first_algorithm.pruned_reactions, first_algorithm.reactions, INPUT_PATH_TXT, OUTPUT_PATH_TXT)
    

    
def run_script(INPUT_PATH_GRAPHML, INPUT_PATH_TXT, OUTPUT_PATH, name="output"):
    # PREPARING THE ALGORITHM
    molecules, reactions, _ = gp.data_parser(INPUT_PATH_GRAPHML)
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
    save_result_graphml(OUTPUT_PATH, first_algorithm, name="testing")
    #OUTPUT_PATH_GRAPHML = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".graphml"
    #gp.to_graphml(INPUT_PATH_GRAPHML, OUTPUT_PATH_GRAPHML, pruned_reactions)
    #
    ## SAVING TO .TXT FILE
    save_result_txt(INPUT_PATH_TXT, OUTPUT_PATH, first_algorithm, name="testing")
    #OUTPUT_PATH_TXT = f"" + OUTPUT_PATH + "/" + name + "_" + first_algorithm.mode + ".txt"
    #ax.run_the_script(pruned_reactions, first_algorithm.reactions, INPUT_PATH_TXT, OUTPUT_PATH_TXT)


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