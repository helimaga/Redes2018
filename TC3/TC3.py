#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TRABAJO COMPUTACIONAL 3
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import community as cm
# https://stackoverflow.com/questions/22070196/community-detection-in-networkx
#import fastcommunity as fg # no lo pude bajar porque la pagina estaba caida
import igraph
# import communityLayout as cl
#import itertools
#import collections
#from random import sample
#import scipy as sp
#from sklearn.linear_model import LinearRegression


pathHeli = '/home/heli/Documents/Redes/Practicas/TC_03/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/tc03/'
pathSanti = '/home/santiago/Documentos/RC/tc03/'
pathDocente = '?'

path = pathJuancho

plt.close('all')
plt.rc('text', usetex=False)
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

#%% Metodología Louvain networkx

dol_part_Louvain = cm.best_partition(dolphins)

#%% Metodología Infomap
#https://www.youtube.com/watch?v=mO0J_H4YLJA
from igraph import *
dolphinsI = Graph.Read_GML(path + 'dolphins.gml')
dol_part_IMap0 = list(dolphinsI.community_infomap())


dol_part_IMap = {}
vs = VertexSeq(dolphinsI)
for n in range(len(dolphins.nodes())):
    vStr = vs[n]['label']
    for comm in range(len(dol_part_IMap0)):
        if n in dol_part_IMap0[comm]:
            dol_part_IMap[vStr] = comm


#%% Metodología Fast Greedy (max Modularity) networkx
# de: https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html

# para importar modulos provenientes de la version trial de networkx...
import sys
sys.path.append(path)
import modularity_max

dol_part_FGreedy0 = list(modularity_max.greedy_modularity_communities(dolphins))
dol_part_FGreedy0 = [list(x) for x in dol_part_FGreedy0]

dol_part_FGreedy = {}
for n in list(dolphins.nodes()):
    for comm in range(len(dol_part_FGreedy0)):
        if n in dol_part_FGreedy0[comm]:
            dol_part_FGreedy[n] = comm
            
#%% Metodología Newman-Girvan (Edge Betweenness) networkx
# http://materias.df.uba.ar/redesa2018c2/files/2018/10/15_Clusters_2.pdf... 
# = Ejercicio TC1.2c  ... pero continuar la division un paso mas? (4 grupos?
# Hay método implementado?

dol_part_NewGir0 = list(nx.algorithms.community.centrality.girvan_newman(dolphins))

dol_part_NewGirAll = [] #lista de diccionarios, cada diccionario es una particion
for numPart in range(len(dol_part_NewGir0)): # len(dol_part_NewGir0) = N-1
    dol_part_NewGir = {}
    for n in list(dolphins.nodes()):
        for comm in range(len(dol_part_NewGir0[numPart])):
            if n in dol_part_NewGir0[numPart][comm]:
                dol_part_NewGir[n] = comm
    dol_part_NewGirAll.append(dol_part_NewGir)
#%% SilhouetteJuancho
def silhouetteJuancho(graph,commPartition,outputOpt):
    numComm = max(commPartition.values())+1
    silhouette = []
    silhouetteAvg = []
    for nSource in graph.nodes():
        cnk = []
        for comm in range(numComm):
            nNodesComm = 0
            dST = 0
            for nTarget in graph.nodes():
                if (commPartition[nTarget] == comm) & (nTarget != nSource):
                    nNodesComm += 1
                    dST += nx.shortest_path_length(graph,source=nSource, target=nTarget)
            nNodesComm = max(nNodesComm,1) # por si nSource es solitario
            dST = max(dST,1) # por si nSource es solitario
            cnk.append(dST/nNodesComm)
        commSource = commPartition[nSource]
        an = cnk[commSource]
        del cnk[commSource]
        bn = min(cnk)
        s0 = (bn-an)/max(an,bn)
        silhouette.append(s0)
    silhouetteAvg = np.mean(silhouette)
    if outputOpt == 'all':
        return silhouette
    elif outputOpt == 'mean':
        return silhouetteAvg
#%% SilhouetteJuancho: Testeando con grafo simple
G = nx.Graph()
G.add_nodes_from([1,2,3,4,5,6])
G.add_edges_from([(1,2),(2,3),(3,1),(1,4),(4,5),(5,6),(6,4)])

GPart = {}
GPart[1] = 0
GPart[2] = 0
GPart[3] = 0
GPart[4] = 1
GPart[5] = 1
GPart[6] = 1

pos={}
pos[1] = np.array((0,1))
pos[2] = np.array((1,2))
pos[3] = np.array((-1,2))
pos[4] = np.array((0,-1))
pos[5] = np.array((1,-2))
pos[6] = np.array((-1,-2))


colour = ['r', 'y', 'b', 'c', 'm']
G_part_col = {}

for n in GPart: 
    G_part_col[n] = colour[GPart[n]]

G_part_col = list(G_part_col.values())

plt.figure()
nx.draw(G,
        pos,
        width=1,
        edge_color = 'm',
        node_color= G_part_col, 
        node_size=70,
        font_size=20,
        with_labels=True,
       )

silTest = str(silhouetteJuancho(G,GPart,'all'))
silMeanTest = str(silhouetteJuancho(G,GPart,'mean'))

plt.suptitle('Testeando Silueta')
plt.text(0, 0.6, 'Sihouette = ' + silTest, fontsize=12)
plt.text(0, 0.5, 'Sihouette Mean = ' + silMeanTest, fontsize=12)
plt.show()

#%% Modularity & Silhouette

silhouetteLouvain = silhouetteJuancho(dolphins,dol_part_Louvain,'mean')
silhouetteFGreedy = silhouetteJuancho(dolphins,dol_part_FGreedy,'mean')
silhouetteIMap = silhouetteJuancho(dolphins,dol_part_IMap,'mean')

modLouvain = cm.modularity(dol_part_Louvain,dolphins)
modFGreedy = cm.modularity(dol_part_FGreedy,dolphins)
modIMap = cm.modularity(dol_part_IMap,dolphins)

silhouetteNewGirAll = []
modNewGirAll = []

for numPart in range(len(dol_part_NewGirAll)):
    silhouetteNewGirPart = silhouetteJuancho(dolphins,dol_part_NewGirAll[numPart],'mean')
    silhouetteNewGirAll.append(silhouetteNewGirPart)
    
    modNewGirPart = cm.modularity(dol_part_NewGirAll[numPart],dolphins)
    modNewGirAll.append(modNewGirPart)

#%% Modularity & Silhouette. Intercomp

plt.subplot(211)

plt.scatter(max(dol_part_Louvain.values())+1,modLouvain, color='g',label='Louvain')
plt.axhline(y=modLouvain, color='g',LineWidth=0.5)
plt.scatter(max(dol_part_FGreedy.values())+1,modFGreedy, color='b',label='Fast Greedy')
plt.axhline(y=modFGreedy, color='b',LineWidth=0.5)
plt.scatter(max(dol_part_IMap.values())+1,modIMap, color='y',label='Info-Map')
plt.axhline(y=modIMap, color='y',LineWidth=0.5)
plt.plot(range(2,len(modNewGirAll)+2),np.array(modNewGirAll),'.-r',label='Newman-Girvan')
plt.xticks(range(0,65))
plt.grid(True,color='k',axis='x',linestyle='--', linewidth=0.5)
plt.xlabel('Number of communities')
plt.ylabel('Modularity')
plt.legend()

plt.subplot(212)

plt.scatter(max(dol_part_Louvain.values())+1,silhouetteLouvain, color='g')
plt.axhline(y=silhouetteLouvain, color='g',LineWidth=0.5)
plt.scatter(max(dol_part_FGreedy.values())+1,silhouetteFGreedy, color='b')
plt.axhline(y=silhouetteFGreedy, color='b',LineWidth=0.5)
plt.scatter(max(dol_part_IMap.values())+1,silhouetteIMap, color='y')
plt.axhline(y=silhouetteIMap, color='y',LineWidth=0.5)
plt.plot(range(2,len(modNewGirAll)+2),np.array(silhouetteNewGirAll),'.-r')
plt.xticks(range(0,65))
plt.grid(True,color='k',axis='x',linestyle='--', linewidth=0.5)
plt.xlabel('Number of communities')
plt.ylabel('Silhouette')

plt.show()
#%% Determinacion de Newman-Girvan
dol_part_NewGir = dol_part_NewGirAll[3]

#%%

colourMethod = [['m', 'g', 'k', 'b', 'r','c'],
                ['r', 'g', 'b', 'k', 'm','c'],
                ['g', 'r', 'b', 'm', 'k','c'],
                ['m', 'g', 'r', 'b', 'c','k']]


methods = [dol_part_Louvain,dol_part_FGreedy,dol_part_IMap,dol_part_NewGir]
methodsStr = ['Louvain','Fast Greedy','Info-Map','Newman-Girvan']

for m in range(len(methods)):
    colour = colourMethod[m]
    dol_part_col = {} 
    for n in methods[m]: 
        dol_part_col[n] = colour[methods[m][n]]

    dol_part_col = list(dol_part_col.values())
    pos = nx.kamada_kawai_layout(dolphins)
    
#    pos=cl.community_layout(dolphins,methods[m])
#        
#    pos=cl.community_layout(dolphins,dol_part_Louvain)
    
    ax = plt.subplot(141+m)
    nx.draw(dolphins,
            pos,
            width=0.1,
            edge_color = 'm',
            node_color= dol_part_col, 
            node_size=50,
            font_size=10,
            with_labels=False,
           )
    ax.set_title(methodsStr[m])
    
plt.show()
#%% k-clique Percolation Method
#%% 1: Encontrar todos los k-clicks de delfines, guardarlos en el dfClick
import networkx.algorithms.clique as click

clicks = list(click.find_cliques(dolphins))
dfClick = pd.DataFrame()
# df of k-clicks
k = 4
c = -1
nodeStr = []
for nck in range(len(clicks)):
    if len(clicks[nck])==k:
        c+=1
        for n in range(k):
            if c==0:
                nodeStr.append('Node' + str(n))
            dfClick.loc[str(c),nodeStr[n]] = clicks[nck][n]
        dfClick.loc[str(c),'Community'] = int(0)
NClicks = c+1
print(dfClick)
#%% 2: Construir un grafo donde cada nodo es un k-click y donde existe enclace
# entre dos clicks si poseen al menos k-1 nodos en comun.
g0 = 0
clickEdges = []
for p in range(NClicks):
    for q in range(p+1,NClicks):
        p0 = list(dfClick.loc[str(p),nodeStr])
        q0 = list(dfClick.loc[str(q),nodeStr])
        r = p0 + q0
        commnodes = len(r)-len(set(r))
        if commnodes>=k-1:
            clickEdges.append([p,q])
    
kclickG = nx.Graph()
kclickG.add_nodes_from(list(range(NClicks)))
kclickG.add_edges_from(clickEdges)

nx.draw(kclickG)

kClickCC = list(nx.connected_component_subgraphs(kclickG))

nodesComm = []
for comm in range(len(kClickCC)):
    clicksComm = list(kClickCC[comm].nodes())
    nodesComm0 = []
    for n in clicksComm:
        nodesComm0 += list(dfClick.loc[str(n),nodeStr])
    nodesComm0 = list(set(nodesComm0))
    nodesComm.append(nodesComm0)

#%%
# plt.sca(axs[0])
#plt.figure()
#nx.draw(dolphins,
#        width=1,
#        edge_color = 'm',
#        node_color= 'k', 
#        node_size=1,
#        font_size=20,
#        with_labels=True,
#       )
#
#plt.suptitle('Red Delfines')
#plt.show()
#
#
#
#
#dol_part_col = {}
#
#for n in dol_part_Louvain: 
#    dol_part_col[n] = colour[dol_part_Louvain[n]]
#
#dol_part_col = list(dol_part_col.values())
#
#plt.sca(axs[1])
##plt.figure()
#nx.draw(dolphins,
#        pos,
#        width=1,
#        edge_color = 'c',
#        node_color= dol_part_col, 
#        node_size=200,
#        font_size=20,
#        with_labels=False,
#       )
#
#plt.suptitle('Red Delfines')
#plt.show()
#%% Metodología Fast Greedy networkx, no sabemos que hace

#dol_part_greedy = nx.algorithms.tree.branchings.greedy_branching(dolphins)
#plt.figure()
#nx.draw(dol_part_greedy,
#        width=1,
#        edge_color = 'c',
#        node_color= dol_part_col, 
#        node_size=200,
#        font_size=20,
#        with_labels=False,
#       )
#%% Metodología Fast Greedy igraph... no funciona elgraficado...
#
#dolphinsI = igraph.Graph.TupleList(dolphins.edges(), directed=True)
#dolphinsI.to_undirected()
#dol_part_greedy = dolphinsI.community_fastgreedy()
#print(dol_part_greedy)
#
#
##color = list(np.random.choice(range(256), size=3))
##print(color)
#
#i = dolphinsI.community_infomap()
#colors = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00","#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"]
#dolphinsI.vs['color'] = [None]
#for clid, cluster in enumerate(i):
#    for member in cluster:
#        dolphinsI.vs[member]['color'] = colors[clid]
#        print(clid)
#dolphinsI.vs['frame_width'] = 0
#igraph.plot(dolphinsI)

# VER:

# https://github.com/taynaud/python-louvain
# https://arxiv.org/pdf/0803.0476.pdf
# https://arxiv.org/pdf/0707.0609.pdf




#%% Metodología infomap
# chequear esto: http://www.mapequation.org/code.html#Linux
#%%

#import networkx as nx
#import matplotlib.pylab as plt
#%matplotlib inline
#
## from rpy2.robjects.packages import importr
## igraph = importr('igraph')
## import pandas as pd
## from rpy2.robjects import r, pandas2ri
## a = pandas2ri.py2ri(nx.to_pandas_adjacency(nxG))
#import igraph
#import os
#import numpy as np
#import rpy2.robjects as robjects
#
##%%
#
#def community(nxG, algorithm, fig_name = "G"):
#    """
#    In:
#        nxG: grafo de networkx.
#        algorithm: string, entre las siguientes opciones: 
#            fast_greedy
#            edge_betweenness
#            louvain
#            infomap
#        fig_name: nombre de la figura que se genera al clsuterizar. Le agrega automaticamente el nombre del algoritmo usado y el nombre del grafo si lo tuviere
#    Out:
#        labels: numpy array con la pertenencia de cada nodo al cluster.
#    
#    """
#    gml_file_name = "G.gml"
#    fig_name += "_"+nxG.name+"_"+algorithm+".svg"
#    nx.write_gml(nxG, gml_file_name)
#    
#    igG = robjects.r('''
#        f <- function(file, algorithm, fig_name){
#            require("igraph")     
#            
#            G <- read_graph(file, "gml")
#            #format = c("edgelist", "pajek", "ncol", "lgl", "graphml","dimacs", "graphdb", "gml", "dl"), ...)
#            
#            if(algorithm == "fast_greedy"){
#                c <- cluster_fast_greedy(G, 
#                    merges = TRUE, 
#                    modularity = TRUE, 
#                    membership = TRUE)
#            }
#            
#            if(algorithm == "edge_betweenness"){
#                c <- cluster_edge_betweenness(G,directed = FALSE,edge.betweenness = TRUE)
#            }
#            
#            if(algorithm == "louvain"){
#                c <- cluster_louvain(G)
#            }
#            
#            if(algorithm == "infomap"){
#                c <- cluster_infomap(G)
#            }
#            
#            svg(fig_name)
#            plot(c, G)
#            dev.off()
#            
#            return(membership(c))
#        }
#    ''')
#    
#    labels = igG(gml_file_name, algorithm, fig_name)
#    os.remove(gml_file_name)
#    return np.array(labels)
#
##%%
#nxG = nx.gnp_random_graph(100, 0.02, seed=None, directed=False)
#nxG.name = "Random"
#
#nx.draw_networkx(nxG)
#plt.axis("off")
#plt.title(nxG.name)
#
##%%
#robjects.r['options'](warn=-1)
#labels_infomap = community(nxG, "infomap")
#labels_infomap = community(nxG, "fast_greedy")
#labels_infomap = community(nxG, "edge_betweenness")
#labels_infomap = community(nxG, "louvain")    

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