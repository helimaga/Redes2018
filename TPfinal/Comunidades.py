#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 11 13:44:26 2018

@author: santiago

Descripción:    Archivo de comunidades, basado en las funciones de Juan para el TC3.
"""

#%%
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import community as cm
import sys


#%%

#paths

pathHeli = '/home/heli/Documents/Redes/Practicas/TPs/Redes2018/TPfinal/'
#pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/tc03/'
pathSanti = '/home/santiago/Documentos/RC/Redes2018/TPfinal/'
pathDocente = '?'

#%%
path = pathSanti

sys.path.append(path)
import modularity_max

#%% Particionamos segun cada uno de los metodos (encapsulamos en funcion...)
def Comunidades(red,metodo):
    if metodo=='l':
        # Metodología Louvain networkx
        dol_part= cm.best_partition(red)

    elif metodo=='fg':
    # Metodología Fast Greedy (max Modularity) networkx
    # de: https://networkx.github.io/documentation/latest/reference/algorithms/generated/networkx.algorithms.community.modularity_max.greedy_modularity_communities.html
    
    # para importar modulos provenientes de la version trial de networkx...
        dol_part_FGreedy0 = list(modularity_max.greedy_modularity_communities(red))
        dol_part_FGreedy0 = [list(x) for x in dol_part_FGreedy0]
        
        dol_part = {}
        for n in list(red.nodes()):
            for comm in range(len(dol_part_FGreedy0)):
                if n in dol_part_FGreedy0[comm]:
                    dol_part[n] = comm
    elif metodo=='im':
        # Metodología Infomap
        #https://www.youtube.com/watch?v=mO0J_H4YLJA
        redI = Graph.Read_GML(path + 'dolphins.gml') # VER ESTE PATH.
        red_part_IMap0 = list(redI.community_infomap())
        
        red_part = {}
        vs = VertexSeq(redI)
        for n in range(len(red.nodes())):
            vStr = vs[n]['label']
            for comm in range(len(red_part_IMap0)):
                if n in red_part_IMap0[comm]:
                    red_part[vStr] = comm
    elif metodo=='ng':
        # Metodología Newman-Girvan (Edge Betweenness) networkx
        # http://materias.df.uba.ar/redesa2018c2/files/2018/10/15_Clusters_2.pdf... 
        # = Ejercicio TC1.2c  ... pero continuar la division un paso mas? (4 grupos?
        # Hay método implementado?
        
        red_part_NewGir0 = list(nx.algorithms.community.centrality.girvan_newman(red))
        
        red_part = [] #lista de diccionarios, cada diccionario es una particion
        for numPart in range(len(red_part_NewGir0)): # len(red_part_NewGir0) = N-1
            red_part_NewGir = {}
            for n in list(red.nodes()):
                for comm in range(len(red_part_NewGir0[numPart])):
                    if n in red_part_NewGir0[numPart][comm]:
                        red_part_NewGir[n] = comm
            red_part.append(red_part_NewGir)
    return red_part
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