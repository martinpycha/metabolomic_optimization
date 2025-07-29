import xml.etree.ElementTree as ET
import collections
from enum import Enum
import html


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
        self.name = name.strip()
        self.hasSource, self.hasTarget = False, False
        self.is_initial_reactant, self.is_final_product = False, False
        self.isPartOfConnectedGraph = False
    def __str__(self):
        return f"{self.name}"
    def __eq__(self, other):
        if not isinstance(other, molecule):
            return NotImplemented
        else:
            return self.name == other.name
    def __hash__(self):
        return hash(self.name)
    def __repr__(self):
        return self.__str__()
    def setInitBools(self):
        self.hasSource, self.hasTarget = False, False
        
class reaction:
    def __init__(self, id, source, target, equation, value, stderr):
        self.id = id
        self.source = source
        self.target = target
        self.equation = equation.strip()
        # TODO TODO
        # tady trochu zapomenu na přesný čísla, protože ty asi dělaji neplechu - 
        #if (value < 0.001):
        #    self.value = 0
        #else:
        #    self.value = value
        self.value = value
        self.stderr = stderr
        self.isRelevant = True  # at the beginning, all reactions are considered to be relevant
        self.measure_rel = measurement_relevance.UNDECIDED
        self.chem_rel = chemical_relevance.UNDECIDED
        self.stat_rel = statistical_relevance.UNDECIDED
    def __str__(self):
        return self.equation + f"\t\t value: {self.value}" + "\n"
        #first = True
        #reaction_string = f""
        #for reactant in self.reactants:
        #    if first:
        #        reaction_string += f"{str(reactant)}"
        #        first = False
        #    else: 
        #        reaction_string += f" + {str(reactant)}"
        #reaction_string += f" -> "
        #first = True
        #for product in self.products:
        #    if first:
        #        reaction_string += f"{str(product)}"
        #        first = False
        #    else: 
        #        reaction_string += f" + {str(product)}"
        ## With or without value?
        #reaction_string += f"\t\t value = {self.value}"
        ##reaction_string += f", stat_rel: {self.stat_rel.value}"
        #reaction_string += '\n'
        ##reaction_string += f", measure_rel = {self.measure_rel}"
        #return reaction_string

    def __repr__(self):
        return self.__str__()
    
    def __eq__(self, other):
        if self.id == other.id:
            return True
        else:
            return False
    def __hash__(self):
        return hash(str(self.id))
        #return hash(self.equation + str(self.source) + str(self.target))




def data_parser(file_path):
    # due to namespace stuff, this line is important:
    ns = {"graphml": "http://graphml.graphdrawing.org/xmlns"}
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    molecules = set()
    node_id_to_mol = {}
    for node in root.findall(".//graphml:node", ns):
        node_id = node.attrib["id"]
        data_elements = {data.attrib['key']: data.text for data in node.findall("graphml:data", ns)}
        node_label = data_elements["d0"]
        node_conf = data_elements["d1"]
        node_includes = data_elements["d2"]
        node_degree = float(data_elements["d3"])
        mol = molecule(node_label)
        node_id_to_mol[node_id] = mol
        molecules.add(mol)
    
    reactions = []
    
    
    for edge in root.findall(".//graphml:edge", ns):
        source_id = edge.attrib.get("source")
        target_id = edge.attrib.get("target")
        reaction_id = edge.attrib.get("id")
        data_elements = {data.attrib['key']: data.text for data in edge.findall("graphml:data", ns)}
        reaction_name = data_elements["d4"]
        reaction_eq = data_elements["d5"]
        reaction_val = float(data_elements["d6"])
        reaction_stdErr = float(data_elements["d7"])
        #reaction_LB = float(data_elements["d8"])
        #
        #if "d12" not in data_elements:
        #    reaction_UB = float(data_elements["d9"])
        #reaction_URL = data_elements["d10"]
        #reaction_path = data_elements["d11"]
        
        source = node_id_to_mol[source_id]
        target = node_id_to_mol[target_id]
        reac = reaction(reaction_id, source, target, reaction_eq, reaction_val, reaction_stdErr)
        reactions.append(reac)
        
        
    print(f"Number of parsed reactions: {len(reactions)}")
        
    return molecules, reactions, node_id_to_mol
    



def to_graphml(input_file_path, output_path, pruned_reactions):
    # Load the GraphML XML
    tree = ET.parse(input_file_path)
    root = tree.getroot()

    # Find the graph element
    graph_elem = root.find("{http://graphml.graphdrawing.org/xmlns}graph")

    # Create a set of IDs to keep
    reaction_ids_to_keep = set(r.id for r in pruned_reactions)

    # Track how many were actually removed
    removed = 0

    # Iterate over all edge elements
    for edge in list(graph_elem.findall("{http://graphml.graphdrawing.org/xmlns}edge")):
        edge_id = edge.get("id")

        # Remove edge if its ID is NOT in the list of reactions to keep
        if edge_id not in reaction_ids_to_keep:
            try:
                graph_elem.remove(edge)
                removed += 1
                #print(f"✅ Removed edge with id={edge_id}")
            except ValueError:
                print(f"error")

    #print(f"Total edges removed: {removed}")

    # Write the pruned GraphML to file
    tree.write(output_path, encoding='utf-8', xml_declaration=True)
                
#file_path = "/Users/martinpycha/Desktop/Job_AV/metabolomic_optimization/threepath.graphml"
#molecules, reactions = data_parser(file_path)

        

            
                
                
                
                
            
                
        
    
    