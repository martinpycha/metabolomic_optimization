import collections
import pandas as pd


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
    def __init__(self, reactants, products, value):
        self.reactants = reactants
        self.products = products
        self.value = value
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
        #reaction_string += f", value = {self.value}"
        return reaction_string
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
   
def simple_trial():
    reactant1 = molecule("Glc")
    reactant2 = molecule("ATP")
    product1 = molecule("GlcP")
    product2 = molecule("ADP")
    #print(product1)
    #print(reactant2)
    reactants = [reactant1, reactant2]
    products = [product1, product2]
    reaction1 = reaction(reactants, products, 100)
    #print(reaction1)
    
    data = pd.read_excel("metabolic_networks/data/Final_data_theoretical_lowR13_withDG.xls")
    print(data.head())
    reactions = []
    for index, row in data.iterrows():
        react, prod = row['Equation'].split(" -> ")
        reactants = [r.strip() for r in react.split(" + ")]
        products = [p.strip() for p in prod.split(" + ")]
        react_val = row['Value']
        single_reaction = reaction(reactants, products, react_val)
        reactions.append(single_reaction)
        
    #for react in reactions:
    #    print(react) 
    return reactions
        
def identify_precs(reactions):
    precursors = set()
    # adding all the reactants into the set of precursors
    for reaction in reactions:
        for reactant in reaction.reactants:
            precursors.add(reactant)
    # removing all the products from the list of precursors
    for reaction in reactions:
        for product in reaction.products:
            precursors.discard(product)
    # removing everything staring with 0* or having .*ATP/ADP
    for precursor in precursors.copy():
        if precursor.startswith("0*") or \
            "ATP" in precursor or \
            "ADP" in precursor:
            precursors.discard(precursor)
    print(precursors)
    
reactions = simple_trial()
identify_precs(reactions)




#def BFS_filter(reactions):
#    queue = []
#    queue.append