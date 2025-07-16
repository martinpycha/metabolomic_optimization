import collections
import pandas as pd
from enum import Enum
import math

# enum according to the relevance of the measured values
class measurement_relevance(Enum):
    UNDECIDED = 'UNDECIDED'
    UNMEASURED = 'UNMEASURED'    # unmeasured
    LOW = 'LOW'                  # the measurements lower than treshold
    HIGH = 'HIGH'                # these measurements are above the treshold
    
class statistical_relevance(Enum):
    UNDECIDED = 'UNDECIDED'
    LOW = 'LOW'                  # the standard deviation is way too high or the interval is too broad
    HIGH = 'HIGH'                # the measurements bear relevance
    
# enum according to the chemical/biological relevance for the system
class chemical_relevance(Enum):
    UNDECIDED = 'UNDECIDED'
    LOW = 'LOW'                 # these reactions/molecules (?) are not likely to bear importance
    HIGH = 'HIGH'               # likely to be important
    CERTAIN = 'CERTAIN'         # certain, must be kept
    
class pruning_mode(Enum):
    MIN_FOREST = 'Minimum spanning forest'


class molecule:
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return f"{self.name}"
    def __eq__(self, other):
        if not isinstance(self, other):
            return NotImplemented
        else:
            return self.name == other.name
        
class reaction:
    def __init__(self, reactants, products, value, stderr):
        self.reactants = reactants
        self.products = products
        self.value = value
        self.stderr = stderr
        self.isRelevant = True  # at the beginning, all reactions are considered to be relevant
        self.measure_rel = measurement_relevance.UNDECIDED
        self.chem_rel = chemical_relevance.UNDECIDED
        self.stat_rel = statistical_relevance.UNDECIDED
    def __str__(self):
        first = True
        reaction_string = f""
        for reactant in self.reactants:
            if first:
                reaction_string += f"{str(reactant)}"
                first = False
            else: 
                reaction_string += f" + {str(reactant)}"
        reaction_string += f" -> "
        first = True
        for product in self.products:
            if first:
                reaction_string += f"{str(product)}"
                first = False
            else: 
                reaction_string += f" + {str(product)}"
        # With or without value?
        reaction_string += f"\t\t value = {self.value}"
        #reaction_string += f", stat_rel: {self.stat_rel.value}"
        reaction_string += '\n'
        #reaction_string += f", measure_rel = {self.measure_rel}"
        return reaction_string

    def __repr__(self):
        return self.__str__()
    
    def __eq__(self, other):
        if collections.Counter(self.reactants) == collections.Counter(other.reactants) and \
            collections.Counter(self.products) == collections.Counter(other.products):
            return True
        else:
            return False
    # is the compound (reactant) part of the reaction?
    def hasReactant(self, reactant):
        for mol in self.reactants:
            if reactant == mol:
                return True
        return False
    # is the compound (product) part of the reaction?
    def hasProduct(self, product):
        for mol in self.products:
            if product == mol:
                return True
        return False
    
# gets panda dataset of reactions, returns List of reaction objects
def get_reactions(pd_dataset):
    reactions = []
    for index, row in pd_dataset.iterrows():
        react, prod = row['Equation'].split(" -> ")
        reactants = [r.strip() for r in react.split(" + ")]
        products = [p.strip() for p in prod.split(" + ")]
        react_val = row['Value']
        reacterr = row['StdErr']
        single_reaction = reaction(reactants, products, react_val, reacterr)
        reactions.append(single_reaction)
    return reactions
    
    
   
#def simple_trial(data_localization):
    
    # trial
    #reactant1 = molecule("Glc")
    #reactant2 = molecule("ATP")
    #product1 = molecule("GlcP")
    #product2 = molecule("ADP")
    #print(product1)
    #print(reactant2)
    #reactants = [reactant1, reactant2]
    #products = [product1, product2]
    #reaction1 = reaction(reactants, products, 100)
    #print(reaction1)
    
    #data = pd.read_excel("metabolic_networks/data/Final_data_theoretical_lowR13_withDG.xls")
    #data = pd.read_excel(data_localization)
    #print(data.head())
    #reactions = []
    #for index, row in data.iterrows():
    #    react, prod = row['Equation'].split(" -> ")
    #    reactants = [r.strip() for r in react.split(" + ")]
    #    products = [p.strip() for p in prod.split(" + ")]
    #    react_val = row['Value']
    #    reacterr = row['StdErr']
    #    single_reaction = reaction(reactants, products, react_val, reacterr)
    #    reactions.append(single_reaction)
    #    
    ##for react in reactions:
    ##    print(react) 
    #return reactions
        
# gets list of reactions
# returns a set reactants and set of products
    

# gets list of reactions
# returns a set of precursors and set of final products AND returns a set of precursors and final products
def identify_precs_and_final_prod(reactions):
    precursors = set()
    final_products = set()
    # adding all the reactants into the set of precursors
    for reaction in reactions:
        for reactant in reaction.reactants:
            precursors.add(reactant)
        for product in reaction.products:
            final_products.add(product)
    # removing all the products from the list of precursors
    all_precs = precursors
    all_prods = final_products
    
    for precursor in precursors.copy():
        if "ATP" in precursor or "ADP" in precursor:
            precursors.discard(precursor)
    for fin_product in final_products.copy():
        if "ATP" in fin_product or "ADP" in fin_product:
            final_products.discard(fin_product)
            
            
    for reaction in reactions:
        for product in reaction.products:
            precursors.discard(product)
        for reactant in reaction.reactants:
            final_products.discard(reactant)
    # removing everything staring with 0* or having .*ATP/ADP 
    #for precursor in precursors.copy():
    #    if precursor.startswith("0*") or \
    #        "ATP" in precursor or \
    #        "ADP" in precursor:
    #        precursors.discard(precursor)
    
    # removing atp, adp etc
    
    
    
    print("PRINTING PRECURSORS:")
    print(precursors)
    print("PRINTING FINAL PRODUCTS:")
    print(final_products)
    return precursors, final_products, all_precs, all_prods
    
    
def calculate_rel_error(value, stderr):
    return math.inf if value == 0 else stderr/value

def calculate_stat_relevance(reactions):
    RELATIVE_ERROR_TRESHOLD = 100.0
    for react in reactions:
        rel_error = calculate_rel_error(react.value, react.stderr)
        if rel_error > RELATIVE_ERROR_TRESHOLD:
            react.stat_rel = statistical_relevance.LOW
        else:
            react.stat_rel = statistical_relevance.HIGH
            
def print_relevant(reactions):
    #print("PRINTING THE STATISTICALY RELEVANT DATA")
    for react in reactions:
        doPrint = True
        for product in react.products:
            if "Sink" in product:
                doPrint = False
        if doPrint:
            if react.stat_rel == statistical_relevance.HIGH:
                print(react)
                
def print_all(reactions):
    print("PRINTING ALL THE REACTIONS")
    for react in reactions:
        print(react)
        
                
class Pruning:
    def __init__(self, reactions, first_reactants, final_products, all_reacts, all_prods, treshold=None, mode='First??'):
        self.reactions = reactions
        self.first_reactants = first_reactants
        self.final_products = final_products
        self.treshold = treshold
        self.mode = mode
        self.all_reactants = all_reacts
        self.all_products = all_prods
        
    def sort(self):
        sorted_reactions = sorted(self.reactions, key=lambda react: react.value)
        print("PRINTING THE SORTED LIST OF REACTIONS")
        print(sorted_reactions)
        return sorted_reactions
        
    # tady dat Kruskala    
    def prune(self):
        pruned_reactions = []
        # postupne prochazet, pridat, kdyz jeste neni
        # + bude treba zajistit, ze cely graf bude spojity!!!
        
        
        
        
        
        ...
       
   
def run_script(data_address):
    #reactions = simple_trial(data_address)
    data = pd.read_excel(data_address)
    reactions = get_reactions(data)
    calculate_stat_relevance(reactions)
    first_reactants, final_products, all_reacts, all_prods = identify_precs_and_final_prod(reactions)
    first_algorithm = Pruning(reactions, first_reactants, final_products, all_reacts, all_prods)
    first_algorithm.sort()
    
    #print_all(reactions)
    
#run_script("/Users/martinpycha/Desktop/Job_AV/metabolic_networks/data/Final_data_theoretical_lowR13_withDG.xlsx")
# tady je puvodni adresa
#run_script("/Users/martinpycha/Desktop/Job_AV/metabolic_networks/data/Final_data_theoretical_lowR13_withDG.xlsx")
# a tady je nova
run_script("/Users/martinpycha/Desktop/Job_AV/metabolomic_optimization/Final_data_theoretical_lowR13_withDG.xlsx")



    #print_relevant(reactions)
#identify_precs(reactions)





#def BFS_filter(reactions):
#    queue = []
#    queue.append