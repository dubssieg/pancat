"Tests for merge and split"
from os import path
from gfagraphs import Graph

graph: Graph = Graph(f"{path.dirname(__file__)}/toy.gfa",
                     "GFA1.1", with_sequence=True)

graph.merge_segments("1", "2")
graph.merge_segments("4", "5")
graph.save_graph(outgraph := f"{path.dirname(__file__)}/toy_1-2_4-5_merge.gfa")


# Le merge est validé !

graph2: Graph = Graph(outgraph,
                      "GFA1.1", with_sequence=True)

graph2.split_segments("1", ['1', '2'], [(0, 3), (3, 5)])
print(graph2.segments)
print(graph2.get_path_list()[0].datas['path'])

graph2.save_graph(f"{path.dirname(__file__)}/toy_1-2_resplit.gfa")
