import networkx as nx
import matplotlib.pyplot as plt

G = nx.read_graphml("metabolic_networks/assets/threepath.graphml")
nx.draw(G, with_labels=True)
plt.show()

