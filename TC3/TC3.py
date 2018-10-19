#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRABAJO COMPUTACIONAL 3
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import itertools
#import collections
#from random import sample
#import scipy as sp
#from sklearn.linear_model import LinearRegression


pathHeli = '/home/heli/Documents/Redes/Practicas/tc03/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/tc03/'
pathSanti = '/home/santiago/Documentos/RC/tc03/'
pathDocente = '?'

path = pathHeli

plt.close('all')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

TitleSize=20
AxisLabelSize=20
LegendSize=15
NumberSize=12
LabelSize=8 # etiquetas de los nodos
NodeSize=50 # tamaño de los nodos
def ldata(archivo):

    f=open(archivo)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)	
    return data

#%%
'''
1) Considere la red social de 62 delfines de Nueva Zelanda (dolphins.txt, dolphins.gml,
dolphinsGender.txt).
'''

dolphins = nx.read_gml(path + 'dolphins.gml') # red
dolphinsGender = np.array(ldata(path+'dolphinsGender.txt')) # numpy array con los generos

# En este loop, identificamos el genero de cada delfin nodo por nodo:
for n in dolphins.nodes():
    dolphins.nodes[n]["gender"] = dolphinsGender[dolphinsGender[:,0]==n,1][0]

# Para corroborar, hacemos un print del cada nodo con su respectivo genero
print(dolphins.nodes("gender"))

#numero de nodos y enlaces de la red dolphins
nodesTot=dolphins.number_of_nodes()
edgesTot=dolphins.number_of_edges()

#%%
'''
a. Encuentre la partición en clusters de esta red utilizando la metodología Louvain, infomap,
fast_greedy y edge_betweenness. Visualice los resultados gráficamente.
'''

# VER:

# https://github.com/taynaud/python-louvain
# https://arxiv.org/pdf/0803.0476.pdf
# https://arxiv.org/pdf/0707.0609.pdf


#%%
'''
b. Caracterice las particiones obtenidas en términos de modularidad y silouhette de cada
partición. Compare con valores esperados en redes recableadas y establezca si tiene derecho a
llamar modular a esta red.
'''

#%%
'''
c. Caracterice cuantitativamente el acuerdo entre las particiones obtenidas utilizando uno o más
de los observables vistos en clase.
'''

#%%
'''
d. Analice cuantitativamente la relación entre el género de los delfines y la estructura de
comunidades del grupo. Puede utilizar para ello, por ejemplo, tests de sobre-representación
y/o sub-representacion. Qué hipótesis puede aventurar sobre propiedades comportamentales
de este grupo de delfines a partir de lo encontrado?
'''

#%%
'''
2) [optativo] Implemente un algoritmo de reconocimiento de comunas basado en la metodología de
percolación de cliques. Qué individuos son los más sociables de la comunidad?
'''