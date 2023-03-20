from gfagraphs import Graph, Segment
from networkx import MultiDiGraph, draw
import matplotlib.pyplot as plt
from rich import print

graph_path: str = "/udd/sidubois/Documents/Code/data/toy_examples/toy.gfa"
graph_version: str = "GFA1.1"
threshold: int = 5

graph_gfagraphs: Graph = Graph(graph_path, graph_version, with_sequence=True)
node_set: set = set()
networkx_graph_before: MultiDiGraph = graph_gfagraphs.compute_networkx()
draw(networkx_graph_before, with_labels=True)
plt.show()

# Display each node length
print('Before :')
print([f"{node.datas['name']} > {node.datas['length']}" for node in graph_gfagraphs.segments])
print('\n'.join(
    [f"{p.datas['name']} : {p.datas['path']}" for p in graph_gfagraphs.get_path_list()]))
print()

# Ne marche que pour les segments où un seul chemin passe par eux, et ne nécessitant pas une cascade de merge
# Si on contrôle le nb de chemins par lequel passe le segment, et que c'est sup. à 1, il faut éviter d'éditer le noeud à la volée.
# Ajouter un champ temporaire au noeud d'info de merge, puis vérifier rien conflictuel ?
for pos, wayline in enumerate(graph_gfagraphs.get_path_list()):
    waypath: list = [wayline.datas['path'][0]]
    for idx, (current_node_name, orientation) in enumerate(wayline.datas['path'][1:], start=1):
        if (current_node := graph_gfagraphs.get_segment(current_node_name)).datas['length'] <= threshold:
            # We need to merge the two nodes
            previous_node: Segment = graph_gfagraphs.get_segment(
                (previous_node_name := wayline.datas['path'][idx-1][0]))
            # Editing the sequence
            current_node.datas['seq'] = previous_node.datas['seq'] + \
                current_node.datas['seq']
            current_node.datas['length'] = len(current_node.datas['seq'])
            # We then need to edit path
            waypath: list = waypath[:-1] + [(current_node_name, orientation)]
        else:
            # We do not need to edit, but we must complete the path
            waypath: list = waypath + [(current_node_name, orientation)]
    wayline.datas['path'] = waypath
    node_set.update(set([path_node for path_node, _ in waypath]))

# In the end, we need to suppress nodes that aren't embed in paths anymore

plt.clf()
# Display each node length
print('After :')
print([f"{node.datas['name']} > {node.datas['length']}" for node in graph_gfagraphs.segments])
print('\n'.join(
    [f"{p.datas['name']} : {p.datas['path']}" for p in graph_gfagraphs.get_path_list()]))
graph_gfagraphs.graph = MultiDiGraph()
networkx_graph_after: MultiDiGraph = graph_gfagraphs.compute_networkx()
draw(networkx_graph_after, with_labels=True)
plt.show()
