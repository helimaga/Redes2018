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
import numpy as np
import scipy
import networkx.algorithms.clique as click
import itertools

from igraph import *
import sys
import modularity_max



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

b. Caracterice las particiones obtenidas en términos de modularidad y silouhette de cada
partición. Compare con valores esperados en redes recableadas y establezca si tiene derecho a
llamar modular a esta red.
'''


#%% Particionamos segun cada uno de los metodos (encapsulamos en funcion...)
def clusteringDolphins(metodo):
    if metodo=='l':
        # Metodología Louvain networkx
        dol_part= cm.best_partition(dolphins)

    elif metodo=='fg':
    # Metodología Fast Greedy (max Modularity) networkx
    # de: https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html
    
    # para importar modulos provenientes de la version trial de networkx...
        sys.path.append(path)
        dol_part_FGreedy0 = list(modularity_max.greedy_modularity_communities(dolphins))
        dol_part_FGreedy0 = [list(x) for x in dol_part_FGreedy0]
        
        dol_part = {}
        for n in list(dolphins.nodes()):
            for comm in range(len(dol_part_FGreedy0)):
                if n in dol_part_FGreedy0[comm]:
                    dol_part[n] = comm
    elif metodo=='im':
        # Metodología Infomap
        #https://www.youtube.com/watch?v=mO0J_H4YLJA
        dolphinsI = Graph.Read_GML(path + 'dolphins.gml')
        dol_part_IMap0 = list(dolphinsI.community_infomap())
        
        dol_part = {}
        vs = VertexSeq(dolphinsI)
        for n in range(len(dolphins.nodes())):
            vStr = vs[n]['label']
            for comm in range(len(dol_part_IMap0)):
                if n in dol_part_IMap0[comm]:
                    dol_part[vStr] = comm
    elif metodo=='ng':
        # Metodología Newman-Girvan (Edge Betweenness) networkx
        # http://materias.df.uba.ar/redesa2018c2/files/2018/10/15_Clusters_2.pdf... 
        # = Ejercicio TC1.2c  ... pero continuar la division un paso mas? (4 grupos?
        # Hay método implementado?
        
        dol_part_NewGir0 = list(nx.algorithms.community.centrality.girvan_newman(dolphins))
        
        dol_part = [] #lista de diccionarios, cada diccionario es una particion
        for numPart in range(len(dol_part_NewGir0)): # len(dol_part_NewGir0) = N-1
            dol_part_NewGir = {}
            for n in list(dolphins.nodes()):
                for comm in range(len(dol_part_NewGir0[numPart])):
                    if n in dol_part_NewGir0[numPart][comm]:
                        dol_part_NewGir[n] = comm
            dol_part.append(dol_part_NewGir)
    return dol_part
#%% Definimos Silhouette
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
#%% Testeamos Silhouette
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

#%% Calculamos Modularidad y silhouette para todas las particiones

dol_part_Louvain = clusteringDolphins('l')
dol_part_FGreedy = clusteringDolphins('fg')
dol_part_IMap = clusteringDolphins('im')
dol_part_NewGirAll = clusteringDolphins('ng')

silhouetteLouvain = silhouetteJuancho(dolphins,dol_part_Louvain,'mean')
silhouetteFGreedy = silhouetteJuancho(dolphins,dol_part_FGreedy,'mean')
silhouetteIMap = silhouetteJuancho(dolphins,dol_part_IMap,'mean')

modLouvain = cm.modularity(dol_part_Louvain,dolphins)
modFGreedy = cm.modularity(dol_part_FGreedy,dolphins)
modIMap = cm.modularity(dol_part_IMap,dolphins)

silhouetteNewGirAll = []
modNewGirAll = []

for numPart in range(len(clusteringDolphins('ng'))):
    silhouetteNewGirPart = silhouetteJuancho(dolphins,dol_part_NewGirAll[numPart],'mean')
    silhouetteNewGirAll.append(silhouetteNewGirPart)
    
    modNewGirPart = cm.modularity(dol_part_NewGirAll[numPart],dolphins)
    modNewGirAll.append(modNewGirPart)

#%% Intercomparamos graficamente modularidad y silhouette

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
#%% Se desprende que la mejor particion Newman-Girvan es la que considera 4 clusters
dol_part_NewGir = dol_part_NewGirAll[3]
#%% Graficamos las particiones de los cuatro metodos (NG = 4 clusters)...

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
#%% Calculamos modularidad sobre redes recableadas a la Maslov
modLouvainS = []
modFGreedyS = []
modIMapS = []
modNewGirS = []

modNewGir = modNewGirAll[3]

N = 500

for shuff in range(0,N):
    print(shuff)
    dolphinsShuff = dolphins.copy()

    dolphinsShuff = nx.double_edge_swap(dolphinsShuff, nswap=90, max_tries=500)
    
    modS = cm.modularity(clusteringDolphins('l'),dolphinsShuff)
    modLouvainS.append(modS)
    
    modS = cm.modularity(clusteringDolphins('fg'),dolphinsShuff)
    modFGreedyS.append(modS)
    
    modS = cm.modularity(clusteringDolphins('im'),dolphinsShuff)
    modIMapS.append(modS)
    
    modS = cm.modularity(clusteringDolphins('ng')[3],dolphinsShuff)
    modNewGirS.append(modS)

methodStr = ['Louvain','FGreedy','IMap','NewGir']
#%% Graficamos... se desprende que la red delfines es altamente modular en 
# comparacion con otras redes recableadas aleatoriamente a la Maslov
plt.figure()
for sp in range(len(methodStr)):
    ax = plt.subplot(221 + sp)
    plt.hist(eval('mod' + methodStr[sp] + 'S'), bins=int(np.floor(N/10)), color='g', edgecolor='k', label='Random Maslov Shuffling')
    plt.axvline(eval('mod' + methodStr[sp]), color='k', linestyle='dashed', linewidth=1, label='Original Network')
    plt.legend()
    plt.xlabel('Modularity')
    plt.ylabel('Occurence')
    ax.set_title(methodStr[sp])
plt.show()


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

#%% k-clique Percolation Method
'''
2) [optativo] Implemente un algoritmo de reconocimiento de comunas basado en la metodología de
percolación de cliques. Qué individuos son los más sociables de la comunidad?
'''
#%% Definiciones para graficar
def ciclo(x,y):
    from functools import reduce
    import operator
    import math
    coords = []
    for t in range(len(x)):
        coords.append([x[t],y[t]])
    center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), coords), [len(coords)] * 2))
    ciclo0 = sorted(coords, key=lambda coord: (-135 - math.degrees(math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360)
    xCiclo = []
    yCiclo = []
    for n in range(len(x)):
        xCiclo.append(ciclo0[n][0])
        yCiclo.append(ciclo0[n][1])
    xCiclo.append(ciclo0[0][0])
    yCiclo.append(ciclo0[0][1])        
    cicloF = {'x': np.array(xCiclo), 'y': np.array(yCiclo)}
    return cicloF

pos = nx.kamada_kawai_layout(dolphins)

colour =  plt.rcParams['axes.prop_cycle'].by_key()['color']

#%%
markerComm = ['*','.','v','^','3','4','8','+','X','H']
ax = plt.subplot(221)
nx.draw(dolphins,
        pos,
        width=0.1,
        edge_color = 'k',
        node_color= 'k', 
        node_size=50,
        font_size=10,
        with_labels=False,
       )
ax.set_title('Dolphins')

ksp = 0
for k in range(3,6):
    ksp += 1
    # Encontrar todos los k-clicks de delfines, guardarlos en el dfClick
    clicks = list(click.find_cliques(dolphins))
    dfClick = pd.DataFrame()
    # df of k-clicks
    nodeStr = []
    for v in range(k):
        nodeStr.append('Node' + str(v))
    
    c = -1
    for nck0 in range(len(clicks)):
        if len(clicks[nck0])>=k:
            allKClicks = list(itertools.combinations(clicks[nck0], k))
            allKClicks = [list(elem) for elem in allKClicks]
            for akc in range(len(allKClicks)):
                c+=1
                n0 = -1
                for n in allKClicks[akc]:
                    n0 += 1
                    dfClick.loc[str(c),nodeStr[n0]] = n
    
    NClicks = c+1
    
    # Construir un grafo donde cada nodo es un k-click y donde existe enclace
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
    
    
    kClickCC = list(nx.connected_component_subgraphs(kclickG))
    
    nodesComm = []
    for comm in range(len(kClickCC)):
        clicksComm = list(kClickCC[comm].nodes())
        nodesComm0 = []
        for n in clicksComm:
            nodesComm0 += list(dfClick.loc[str(n),nodeStr])
        nodesComm0 = list(set(nodesComm0))
        nodesComm.append(nodesComm0)
    # Graficar
    ax = plt.subplot(221+ksp)
    for comm in range(len(nodesComm)):
        x = []
        y = []
        for n in nodesComm[comm]:
            x.append(pos[n][0])
            y.append(pos[n][1])
        
        ciclo0 = ciclo(x,y)
        
        x = ciclo0['x']
        y = ciclo0['y']
        
    #    tck, u = scipy.interpolate.splprep([x,y], k = min(len(x)-1,5))
    #    unew = np.arange(0, 1.001, 0.001)
    #    out = scipy.interpolate.splev(unew, tck)
        plt.plot(x,y,markerComm[comm],color=colour[comm],markersize=15)
        plt.plot(x,y,'-',color=colour[comm],linewidth=1)
    nx.draw(dolphins,
            pos,
            width=0.1,
            edge_color = 'k',
            node_color= 'k', 
            node_size=50,
            font_size=10,
            with_labels=False,
           )
    ax.set_title(str(k) + '-cilque Percolation')
plt.show()