#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 10:07:18 2018

@author: gossn
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 14:58:22 2018

@author: heli
"""
#%% Importar modulos

import pandas as pd
import itertools as itr
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62 
from Bio.SubsMat.MatrixInfo import blosum90
import matplotlib.patches as mpatches
import csv
import sys
from pathlib import Path
import pickle
from sklearn import cluster
from sklearn import metrics
import ast
import community
from matplotlib import colors as mcolors
from numpy import random
import scipy.stats as stats
#%% Determinar path e importar modulo casero "Comunidades" (simil TC3)

pathHeli = '/home/heli/Documents/Redes/Practicas/TPs/Redes2018/TPfinal/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/redes2018/TPfinal/'
pathSanti = '/home/santiago/Documentos/RC/Redes2018/TPfinal/'
pathDocente = '?'


path = pathJuancho

sys.path.append(path)
import Comunidades as COM


#%% OPCIONES DE ACHICAMIENTO DE RED

nSample = 50 # si -1, no quita ninguno...
virusSelect = [0,1,2,3,4,5]

#%% CARGAR DATOS VDJDV

virusStr = ['YFV', 'CMV', 'EBV', 'HIV-1', 'InfAV', 'HCV']
virusStr = list( virusStr[i] for i in virusSelect)


TRBVirDf = {} # Diccionario de DataFrames

for i in virusStr:

    # Cargo las tablas con los TCRs asociados a epitopes asociados a cada especie de virus en la lista
    df0 = pd.read_csv(path + 'Datos/'+ i, sep='\t')

    #filtro por las secuencias con Score mayor a 0 (con algun grado de confianza para su anotacion), Gene igual a TRB == T-cell Receptor Beta, Species igual a Homo Sapiens y MHC_class (clase del Complejo Mayor de Histocompatibilidad) igual a I. 
    df1 = df0[(df0.Score > 0) &(df0.Gene == 'TRB') & (df0.Species == 'HomoSapiens') & (df0.MHC_class == 'MHCI')] # 

    #me quedo con las columnas CDR3 (Complementary Determining Region 3) que corresponde a la secuencia de aminoacidos del TRB, Epitope que corresponde a la secuencia del epitope y Epitope_gene que corresponde al gen del cual proviene dicho epitope.
    df2 = df1 [['CDR3', 'Epitope', 'Epitope_gene', 'Epitope_species']]
    if nSample!=-1:
            df2 = df2.sample(n=nSample)
    print(len(df2))

#    # Elimino secuencias CDR3 duplicadas
#    df3 = df2.drop_duplicates(subset = 'CDR3')

    TRBVirDf[i] = df2
#%% Armar un unico DataFrame con todxs lxs virus :P

TRBdf0 = []

for i in virusStr:
    TRBdf0.append(TRBVirDf[i]) # Lista de dataframes

# Dataframe
TRBdf = pd.concat(TRBdf0)

# Reseteo de indices, nombrarlos de 0 a N...
TRBdf = TRBdf.reset_index(drop=True)


#%% Cargar METRICAS

# OBS: Metricas Blosum62 y Blosum90 ya importadas (ver primer bloque)

##### PMHC

# Cargo la matriz de similaridad del paper de Kim et al. (2009) (novelmatrixTCR.pdf) 
peptideMHCmatrix = pd.read_csv(path + 'newmatrix.txt', sep=' ')
peptideMHCmatrix.columns = (list(peptideMHCmatrix)[1:21]+[0])
peptideMHCmatrix = peptideMHCmatrix.iloc[:,0:20]

#para reescalear todos los valores de la matriz (dataframe)
#peptideMHCmatrix = peptideMHCmatrix.multiply(1000).astype(int)

pMHCmatrix = {}

pairsmatrix = list(itr.combinations_with_replacement(list(peptideMHCmatrix),2))
for (i,j) in pairsmatrix:
    pMHCmatrix[(i,j)] = peptideMHCmatrix[i][j]


# OJOOOO!!!!! REEMPLAZAR por la metrica adecuada!!!
from Bio.SubsMat.MatrixInfo import blosum50


##### Inserte su nueva metrica AQUI: 

#%% ELIMINO DATAFRAMES AL PEDO... ASI al apretar en la columna de "Type" en el 
# variable explorer los tienen todos arriba y pipí cucú
del(df0,df1,df2,TRBVirDf,peptideMHCmatrix)
#%% DEFINIR REDES

##### OJOOOO: Reemplazar blosum50 con la metrica correcta

columns = ['nodeType', 'matrixName', 'matrix', 'gop', 'gep', 'computOpt']

dfNetworks = pd.DataFrame([\
                           ['CDR3','blosum62',blosum62,-10,-1,'pairwiseAlign'],\
                           ['CDR3','blosum50',blosum50,-10,-1,'pairwiseAlign'],\
                           ['Epitope','blosum90',blosum90,-10,-1,'pairwiseAlign'],\
                           ['Epitope','pMHCmatrix',pMHCmatrix,-1,-0.1,'pairwiseAlign'],\
                           ], columns=columns)

nNets = dfNetworks.index

#%% COMPUTAR SCORES!
#tomo todos los pares unicos de secuencias y hago pairwise alignment de las 
#mismas para obtener una score

nodesNets = {}
scoresNets = {}

for net in nNets:
    print(net)
    # 1: Elijo los nodos
    nodos = list(TRBdf[dfNetworks.loc[net,'nodeType']])
    # 2: Elimino nodos repetidos
    nodos = list(set(nodos))
    # 3: Los ordeno alfabeticamente
    nodos = sorted(nodos)
    
    nodesNets[net] = nodos
    
    nodePairs = list(itr.combinations(nodos, 2))
    scorelist = []
    
    matrix = dfNetworks.loc[net,'matrix']
    gop = dfNetworks.loc[net,'gop']
    gep = dfNetworks.loc[net,'gep'] 
    computOpt = dfNetworks.loc[net,'computOpt'] 
    
    for (seq1, seq2) in nodePairs:
        if computOpt == 'pairwiseAlign':
            score = pairwise2.align.localds(seq1, seq2, matrix, gop, gep, score_only=True)
        else:
            print('Introduzca su nuevo metodo AQUI')
            score = 1
        scorelist.append(score)
    scoresNets[net] = scorelist

#%% DEFINIR ENLACES A PARTIR DE UMBRAL

perScore = {'CDR3': 0.95, 'Epitope': 0.95}

edgesNets = {}
threshNets = {}
graphsNets = {}


for net in nNets:
    print(net)
    threshNets[net] = np.percentile(scoresNets[net],100*perScore[dfNetworks.loc[net,'nodeType']])
    nodePairs = list(itr.combinations(nodesNets[net], 2))
    selectedPairs = []
    pair = -1
    for (seq1, seq2) in nodePairs:
        pair+=1
        if scoresNets[net][pair] > threshNets[net]:
            selectedPairs.append((seq1, seq2))
    edgesNets[net] = selectedPairs
    
    graphsNets[net] = nx.Graph()
    graphsNets[net].add_edges_from(edgesNets[net])
    graphsNets[net] = max(nx.connected_component_subgraphs(graphsNets[net]), key=len)
    #    graphsNets[net].add_nodes_from(nodesNets[net])
#%% HISTOGRAMAS SCORES + UMBRAL

#hago un histograma con los scores de las secuencias
f = plt.figure()
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()

sp=0
for net in nNets:
    sp+=1
    plt.subplot(220+sp)
    plt.hist(scoresNets[net], bins=500)
    plt.axvline(x=threshNets[net], color='k',LineWidth=0.5)
    #plt.xticks(np.arange(0, 80, step=5))
    plt.xlim((0,1.5*threshNets[net]))
    titStr = dfNetworks.loc[net,'nodeType'] + '; ' + dfNetworks.loc[net,'matrixName']
    plt.title(titStr, fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=10)
    if net>1:
        plt.xlabel('Sequence score', fontsize=20)
    plt.ylabel('Occurrence', fontsize=20)
    plt.show()

#%% CLUSTERIZACION: INFOMAP, LOUVAIN, ...

clusterType = {'InfoMap':{'method':'TC3','acr':'im'},\
               'Louvain':{'method':'TC3','acr':'l'},}

cMethods = list(clusterType.keys())

clusterNets = pd.DataFrame()
clusterNets = clusterNets.astype('object')
idx = -1
for net in nNets:
    for clust in cMethods:
        print(str(net) + ' ' + clust)
        idx+=1
        clusterNets.loc[idx,'IDNet'] = int(net)
        clusterNets.loc[idx,'clusterType'] = clust
        if clusterType[clust]['method'] == 'TC3':
            methodStr = clusterType[clust]['acr']
            clusterNets.loc[idx,'clusterDict'] = [COM.Communities(graphsNets[net],path,methodStr)] # Va como lista, no como diccionario...
        commPartition = clusterNets.loc[idx,'clusterDict'][0]
        clusterNets.loc[idx,'silhouette'] = [COM.silhouetteJuancho(graphsNets[net],commPartition,'all')]
        clusterNets.loc[idx,'silhouetteAvg'] = COM.silhouetteJuancho(graphsNets[net],commPartition,'mean')
        clusterNets.loc[idx,'modularity'] = community.modularity(commPartition,graphsNets[net])
clusterNets['IDNet'] = clusterNets['IDNet'].astype(int)

#%% GRAFOS para InfoMap y Louvain
colores = ['aqua',
 'blueviolet',
 'peru',
 'darkgrey',
 'darkseagreen',
 'violet',
 'red',
 'hotpink',
 'moccasin',
 'green',
 'lime',
 'mediumpurple',
 'sandybrown',
 'orangered',
 'peru',
 'royalblue',
 'skyblue',
 'teal',
 'yellow']

markerComm = ['*','.','v','^','3','4','8','+','X','H']
#colors = list(dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS).keys())
#colores = list(range(len(colors)))
#random.shuffle(colores)

for net in nNets:
    commPart = [\
    list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[0]),'clusterDict'])[0][0],\
    list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[1]),'clusterDict'])[0][0]]
    A,B = COM.dictsValues2Mat(commPart[0],commPart[1])
    commPart[0] = A
    commPart[1] = B
    
    pos = nx.kamada_kawai_layout(graphsNets[net])
        
    plt.figure(net)
    sp = -1
    
    netType = dfNetworks.loc[net,'nodeType']
    metrica = dfNetworks.loc[net,'matrixName']
    
    for clust in cMethods:
        sp+=1
        colNodeList = []
        for n in graphsNets[net].nodes():
            colNodeList.append(colores[commPart[sp][n]])
        ax = plt.subplot(121+sp)
        nx.draw(graphsNets[net],
                pos,
                width=0.2,
                edge_color = 'k',
                node_color= colNodeList, 
                node_size=50,
                font_size=10,
                with_labels=False,
               )
        ax.set_title('Net: ' + netType + '-' + metrica + '\n' + 'Partition: ' + clust)
#%% SILUETAS para InfoMap y Louvain
for net in nNets:
    commPart = [\
    list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[0]),'clusterDict'])[0][0],\
    list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[1]),'clusterDict'])[0][0]]
    A,B = COM.dictsValues2Mat(commPart[0],commPart[1])
    commPart[0] = A
    commPart[1] = B

    plt.figure(net)
    sp = -1
    
    netType = dfNetworks.loc[net,'nodeType']
    metrica = dfNetworks.loc[net,'matrixName']
    for clust in cMethods:
        sp+=1

        dfsillpart = pd.DataFrame()
        silueta0 = list(clusterNets.loc[(clusterNets['IDNet']==net)&(clusterNets['clusterType']==clust),'silhouette'])[0][0]
        silmed = float(clusterNets.loc[(clusterNets['IDNet']==net)&(clusterNets['clusterType']==clust),'silhouetteAvg'])
        print(silmed)
        for n in graphsNets[net].nodes():
            dfsillpart.loc[n,'community'] = commPart[sp][n]
            dfsillpart.loc[n,'color'] = colores[commPart[sp][n]]
            dfsillpart.loc[n,'silhouette'] = silueta0[n]
        dfsillpart = dfsillpart.sort_values(by=['community', 'silhouette'], ascending=[True, False])
        
        silueta1 = list(reversed(dfsillpart['silhouette']))
        color1 = list(reversed(dfsillpart['color']))
        nodos1 = list(reversed(dfsillpart.index))
        
        
        color2 = []
        for n in graphsNets[net].nodes():
            color2.append(dfsillpart.loc[n,'color'])
        
        ax = plt.subplot(121+sp)
        
        plt.barh(range(len(nodos1)), silueta1, align='center',color=color1)
        plt.axvline(x=silmed,color='k',linestyle='--')
        plt.axvline(x=0,color='k',linestyle='-',linewidth=0.5)
        plt.text(silmed+0.01,1,'mean= ' + str(round(silmed*100)/100),fontsize=12,color='k')
        plt.yticks(range(len(nodos1)),'')
        plt.xlabel('silueta')
        plt.xlim(-1,1)
        ax.set_title('Net: ' + netType + '-' + metrica + '\n' + 'Partition: ' + clust)
        plt.show()
        plt.draw()
del(dfsillpart)
#%% TESTS de FISHER 2x2 (SIGUE los IDS del DF clusterNets)

# Pueden elegir el criterio separador que les parezca acá (columnas de TRBdf)
sepCrit = 'Epitope_species'
pThresh = 0.05 # p Muy exigente

contMatrices = []
fisherTests = []


for net in nNets:
    sepDict = {}
    comPart = {}
    N = len(graphsNets[net].nodes())
    for n in graphsNets[net].nodes():
        sepDict[n] = list(TRBdf.loc[TRBdf[dfNetworks.loc[net,'nodeType']]==n,sepCrit])[0]
    comPart = list(clusterNets.loc[(clusterNets['IDNet']==net) & (clusterNets['clusterType']==cMethods[0]),'clusterDict'])[0][0]
    
    sepCritVal = list(sepDict.values())
    comPartVal = list(comPart.values())
    sepCritList = list(set(sepCritVal))
    comPartList = list(set(comPartVal))

    for clust in cMethods:
        contMatrix = pd.DataFrame(columns=sepCritList)
        fisherTest = pd.DataFrame(columns=sepCritList)
        for com in comPartList:
            A = comPartVal.count(com)
            for sepval in sepCritList:
                Nsepval = sepCritVal.count(sepval)
                Asepval = 0
                for n in graphsNets[net].nodes():
                    if (sepDict[n] == sepval) & (comPart[n] == com):
                        Asepval+=1
                contMatrix.loc[com,sepval] = Asepval
                oddsratio, pvalue = stats.fisher_exact([[Asepval,A-Asepval], [Nsepval-Asepval,N-A-(Nsepval-Asepval)]],'less')
                if pvalue > 1-pThresh:
                    fisherTest.loc[com,sepval] = '+' + sepval
                else:
                    fisherTest.loc[com,sepval] = 'Neutro'
        contMatrices.append(contMatrix)
        fisherTests.append(fisherTest)
        print(contMatrix)
        print(fisherTest)
        print('\n\n\n')
        
#%%