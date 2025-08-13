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
    def __init__(self, name, id=0):
        self.name = name.strip()
        self.hasSource, self.hasTarget = False, False
        self.is_initial_reactant, self.is_final_product = False, False
        self.isPartOfConnectedGraph = False
        self.id = id
        self.old_names = set()
        self.old_ids = set()
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
    def addOldNameId(self, another_old_name, another_old_id):
        self.old_names.add(another_old_name)
        self.old_ids.add(another_old_id)
        
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
        mol = molecule(node_label, node_id)
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




import xml.etree.ElementTree as ET

def to_graphml_visualize(input_file_path, output_path, pruned_reactions):
    """
    Write a GraphML that:
      - keeps only reactions in pruned_reactions
      - updates edge source/target to the (possibly unified) molecule ids
      - updates edge equation text (d5) from reaction.equation
      - removes nodes not referenced by any kept edge
      - updates molecule node names (d0) to the molecule.name used in pruned_reactions
      - creates missing nodes if some used molecule id has no node in the input
    """
    NS = "http://graphml.graphdrawing.org/xmlns"
    ns = {"g": NS}

    # Parse
    tree = ET.parse(input_file_path)
    root = tree.getroot()
    graph = root.find("g:graph", ns)
    if graph is None:
        raise ValueError("No <graph> element found in GraphML.")

    # Build quick lookups from pruned_reactions
    # Reaction map by GraphML edge id (string)
    reaction_by_id = {str(r.id): r for r in pruned_reactions}

    # Molecules actually used in pruned reactions (id -> molecule object)
    id_to_mol = {}
    for r in pruned_reactions:
        id_to_mol[str(r.source.id)] = r.source
        id_to_mol[str(r.target.id)] = r.target
    used_mol_ids = set(id_to_mol.keys())

    # Cache existing nodes in the XML by id for quick access
    xml_nodes_by_id = {}
    for node in graph.findall("g:node", ns):
        xml_nodes_by_id[node.get("id")] = node

    # ---------- EDGES ----------
    # Remove edges not in pruned_reactions; update kept edges from reaction objects
    for edge in list(graph.findall("g:edge", ns)):
        eid = edge.get("id")
        react = reaction_by_id.get(eid)
        if react is None:
            # Drop edges not in the pruned result
            graph.remove(edge)
            continue

        # Update source/target to match unified molecule ids
        edge.set("source", str(react.source.id))
        edge.set("target", str(react.target.id))

        # Update equation (assumes d5 is equation field)
        for data_elem in edge.findall("g:data", ns):
            if data_elem.get("key") == "d5":
                data_elem.text = react.equation
                break  # stop after first match

    # ---------- NODES ----------
    # 1) Remove any node not referenced by kept edges (i.e., not in used_mol_ids)
    for node in list(graph.findall("g:node", ns)):
        nid = node.get("id")
        if nid not in used_mol_ids:
            graph.remove(node)
            # also drop from cache so we know it's gone
            xml_nodes_by_id.pop(nid, None)

    # 2) Ensure every used molecule id has a node; if missing, create it
    #    Try to clone from any of its old_ids if present; otherwise create minimal node with d0
    for mid in used_mol_ids:
        if mid in xml_nodes_by_id:
            continue  # node exists already

        mol = id_to_mol[mid]
        new_node = None

        # Try to find an old node to clone (if this molecule has old_ids)
        template = None
        for old_id in getattr(mol, "old_ids", set()):
            template = xml_nodes_by_id.get(str(old_id))
            if template is not None:
                break

        if template is not None:
            new_node = deepcopy(template)
            new_node.set("id", str(mid))
            # update name field (d0)
            for d in new_node.findall("g:data", ns):
                if d.get("key") == "d0":
                    d.text = mol.name
        else:
            # Create a minimal node with just the name
            new_node = ET.Element(f"{{{NS}}}node", {"id": str(mid)})
            d0 = ET.SubElement(new_node, f"{{{NS}}}data", {"key": "d0"})
            d0.text = mol.name

        graph.append(new_node)
        xml_nodes_by_id[str(mid)] = new_node

    # 3) Update names (d0) of all kept/created nodes to match current molecule names
    for mid in used_mol_ids:
        node = xml_nodes_by_id.get(mid)
        if node is None:
            continue
        mol = id_to_mol[mid]
        # set d0 to mol.name
        updated = False
        for d in node.findall("g:data", ns):
            if d.get("key") == "d0":
                d.text = mol.name
                updated = True
                break
        if not updated:
            d0 = ET.SubElement(node, f"{{{NS}}}data", {"key": "d0"})
            d0.text = mol.name

    # Write out
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
                
#file_path = "/Users/martinpycha/Desktop/Job_AV/metabolomic_optimization/threepath.graphml"
#molecules, reactions = data_parser(file_path)

        

            
                
                
                
                
            
                
        
    
    