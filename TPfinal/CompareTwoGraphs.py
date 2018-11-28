"""
Diferencia entre dos grafos.
"""

import networkx as nx
import sys
import matplotlib.pyplot as plt

pathHeli = '/home/heli/Documents/Redes/Practicas/TPs/Redes2018/TPfinal/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/redes2018/TPfinal/'
pathSanti = '/home/santiago/Documentos/RC/Redes2018/TPfinal/'
pathDocente = '?'

#%%
path = pathSanti

plt.close('all')

graph1 = [("b","c"), ("b","d"), ("d","e"), ("e","f"), ("b","f"), ("b","g"), ("f", "g"), ("a","g"), ("a","g"), ("c","d"), ("d", "g"), ("d","h"), ("aa","h"), ("aa","c"), ("f", "h")]
graph2 = [("a","b"), ("b","c"), ("b","d"), ("d","e"), ("e","f"), ("b","f"), ("b","g"), ("f", "g"), ("a","g"), ("a","g"), ("c","d"), ("d", "g"), ("d","h"), ("aa","h"), ("aa","c")]#, ("f", "h")] + [ ("h", "m"), ("m", "l"), ("l", "C"), ("C", "r"), ("a", "k"), ("k", "l"), ("k", "C")] #+ [("z","zz")]

graph1 = nx.Graph(graph1)
graph2 = nx.Graph(graph2)

Grafos=[graph1,graph2]
ColorGrafos=['red','blue']

# Visualizo los dos grafos:
for i in range(len(Grafos)):
    plt.figure()
    nx.draw(Grafos[i],
            width=0.1,
            edge_color = 'k',
            node_color = ColorGrafos[i], 
            node_size=50,
            font_size=30,
            with_labels=True,
            )
    plt.show()

d11=nx.algorithms.similarity.graph_edit_distance(graph1,graph1)
d12=nx.algorithms.similarity.graph_edit_distance(graph1,graph2)

print('Diferencia entre grafos 1 y 1 (i.e. mismo grafo):',d11)
print('Diferencia entre grafos 1 y 2 (i.e. grafos distintos):',d12)

# Hasta donde entendí jugando con los nodos, la distancia aumenta una unidad por cada enlace o nodo en graph2 que no esté en 