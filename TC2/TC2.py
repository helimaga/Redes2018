#TRABAJO COMPUTACIONAL 2#

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
import collections
from random import sample

'''
import scipy as sp
from scipy import stats
from matplotlib import pylab
from sklearn.linear_model import LinearRegression

import matplotlib.patches as mpatches
from statsmodels.stats.proportion import proportions_ztest

import igraph

from matplotlib_venn import venn3, venn3_circles
'''

pathHeli = '/home/heli/Documents/Redes/Practicas/TC_02/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/tc02Data/'
pathSanti = '/home/santiago/Documentos/RC/tc02Data/'
pathDocente = '?'

path = pathSanti

plt.close('all')
plt.rc('text', usetex=False)
plt.rc('font', family='serif')

TitleSize=20
AxisLabelSize=20
LegendSize=11
NumberSize=12
LabelSize=8 # etiquetas de los nodos
NodeSize=50 # tamaño de los nodos

#%%

redesStr = ['Y2H','AP-MS','LIT','LIT_Reguly']
redes = {}

for s in redesStr:
    df = pd.read_csv(path + 'yeast_' + s + '.txt', sep = '\t')
    if s == 'LIT_Reguly':
        net = df[['Bait gene/protein','Hit gene/protein']].values.tolist()
    else:
        net = df.values.tolist()
    redes[s] = nx.Graph(net)

df = pd.read_csv(path + 'Essential_ORFs_paperHe.txt', sep='\t')
esenciales = df['ORF_name'].values.tolist()

#borramos el primer y los dos ultimos elementos de la lista esenciales
del esenciales[0]
del esenciales[len(esenciales)-1]
del esenciales[len(esenciales)-1]

print(esenciales)

#sacamos los espacios al final de cada nombre de la lista
essentials = []
for item in range(len(esenciales)):
    essentials.append(esenciales[item].strip())
    
print(essentials)

#%%

#Caracteristicas de las redes analizadas - Tabla 2 de Zotenko et al.

#primero chequeo que las redes no tengan enlaces repetidos

for s in redesStr:
    if len(list(redes[s].edges()))==len(np.unique(list(redes[s].edges()), axis=0)):
        print("No hay enlaces repetidos para la red " + s)
    else:
        print("Hay enlaces repetidos para la red " + s)

iter=list(itertools.combinations(range(len(redesStr)),2))

def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2)) 

ratios = np.zeros((4, 4), float)
np.fill_diagonal(ratios, 1)

for (j,k) in iter:
    red_j = redes[redesStr[j]].copy()
    red_k = redes[redesStr[k]].copy()
    
    commonNodes = intersection(list(red_j), list(red_k))
    red_j.remove_nodes_from(list(nodoi for nodoi in red_j if nodoi not in commonNodes))
    red_k.remove_nodes_from(list(nodoi for nodoi in red_k if nodoi not in commonNodes))
    
    red_jk = nx.intersection(red_j,red_k)
    commonEdges = red_jk.number_of_edges()
    
    ratios[j,k] = commonEdges/(redes[redesStr[j]].number_of_edges())
    ratios[k,j] = commonEdges/(redes[redesStr[k]].number_of_edges())

dfB2 = pd.DataFrame(ratios, columns=redesStr)
dfB2.index = redesStr

#%% 

#Caracteristicas de las redes analizadas - Figura 1A  de Zotenko et al.

percentiles = np.linspace(0,100,num=101)

degree_sequence = {}
degree_cutoff = {}
ratio_hubs_essentials = {}

for k in redesStr:
    degree_sequence[k] = sorted([d for n, d in redes[k].degree()])
    deg_cutoff = []
    for i in range(len(percentiles)):
        deg_cutoff.append(np.percentile(degree_sequence[k], i))
    degree_cutoff[k] = deg_cutoff

    ratio = []
    for j in range(len(deg_cutoff)):
        count_hubs = 0
        count_hubs_essentials = 0
        for nodoi in list(redes[k].nodes()):
            if redes[k].degree(nodoi)>=deg_cutoff[j]:
                count_hubs+=1
                if nodoi in essentials:                                                                            
                    count_hubs_essentials+=1
        ratio.append(count_hubs_essentials/count_hubs)
    ratio_hubs_essentials[k]= ratio


x = []
for i in range(len(percentiles)):
    x.append(1 - percentiles[i]/100)


plt.figure()
for k in redesStr:
    plt.plot(x, ratio_hubs_essentials[k], linewidth=3, label=k)
plt.xticks(np.arange(0, 1.1, step=0.1))
plt.yticks(np.arange(0, 1.1, step=0.1))
plt.tick_params(axis='both', which='major', labelsize=NumberSize)
plt.xlabel('Hub definition cutoff', fontsize=20)
plt.ylabel('Fraction of essential nodes', fontsize=20)   
plt.title('Figure 1A. Relationship between degree and essentiality in the tested networks', fontsize=30)
plt.legend(loc='upper right', fontsize=20) 
plt.grid(axis='both', color='k', linestyle='dashed', linewidth=2, alpha=0.1)
plt.show()


#%%

#Analisis de vulnerabilidad
#Figura 3 de Zotenko et. al. (2008)


for s in redesStr:
    ######################################################
    # Aleatorio:
    print('\nRANDOM:')
    red=redes[s].copy()
    x_random=[]
    y_random=[]
    N0=red.number_of_nodes()
    print('Cantidad inicial de nodos de la red %s:\t%d'%(s,N0))
    
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component=giant.number_of_nodes()
    while largest_component>1:
        randomNode=np.random.choice(list(red.nodes))   # se elige un nodo al azar
        red.remove_node(randomNode)
        giant = max(nx.connected_component_subgraphs(red), key=len)
        largest_component=giant.number_of_nodes()
        N=red.number_of_nodes()
        
        print('\rNodos restantes:\t%4d'%(N),end='')
        x_random.append(1-N/N0)
        y_random.append(largest_component/N)

    ######################################################
    # Grado:
    print('\nDEGREE:')
    red=redes[s].copy()
    x_degree=[]
    y_degree=[]
    N0=red.number_of_nodes()
    print('Cantidad inicial de nodos de la red %s:\t%d'%(s,N0))
    
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component=giant.number_of_nodes()
    while largest_component>1:
        D=dict(nx.degree_centrality(red))
        maxDegreeNode = max(D, key=D.get) # se toma el nodo de máxima centralidad.
        red.remove_node(maxDegreeNode)
        giant = max(nx.connected_component_subgraphs(red), key=len)
        largest_component=giant.number_of_nodes()
        N=red.number_of_nodes()
        
        print('\rNodos restantes:\t%4d'%(N),end='')
        x_degree.append(1-N/N0)
        y_degree.append(largest_component/N)
    
    ######################################################
    # Autovalor:
    print('\nEIGENVECTOR:')
    red=redes[s].copy()
    x_eigenvector=[]
    y_eigenvector=[]
    N0=red.number_of_nodes()
    print('Cantidad inicial de nodos de la red %s:\t%d'%(s,N0))
    
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component=giant.number_of_nodes()
    while largest_component>1:
        D=dict(nx.eigenvector_centrality(red,tol=1e-03))  # Se relaja la tolerancia (por defecto 1e-06) porque si no, no converge.
        maxDegreeNode = max(D, key=D.get) # se toma el nodo de máxima centralidad.
        red.remove_node(maxDegreeNode)
        giant = max(nx.connected_component_subgraphs(red), key=len)
        largest_component=giant.number_of_nodes()
        N=red.number_of_nodes()
        
        print('\rNodos restantes:\t%4d'%(N),end='')
        x_eigenvector.append(1-N/N0)
        y_eigenvector.append(largest_component/N)
    
    
    ######################################################
    # Subgrafos:
    print('\nSUBGRAPH:')
    red=redes[s].copy()
    x_subgraph=[]
    y_subgraph=[]
    N0=red.number_of_nodes()
    print('Cantidad inicial de nodos de la red %s:\t%d'%(s,N0))
    
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component=giant.number_of_nodes()
    while largest_component>1:
        D=dict(nx.subgraph_centrality(red))
        maxDegreeNode = max(D, key=D.get) # se toma el nodo de máxima centralidad.
        red.remove_node(maxDegreeNode)
        giant = max(nx.connected_component_subgraphs(red), key=len)
        largest_component=giant.number_of_nodes()
        N=red.number_of_nodes()
        
        print('\rNodos restantes:\t%4d'%(N),end='')
        x_subgraph.append(1-N/N0)
        y_subgraph.append(largest_component/N)
    
    ######################################################
    # Camino más corto:
    print('\nSHORTEST PATH:')
    red=redes[s].copy()
    x_shortest_path=[]
    y_shortest_path=[]
    N0=red.number_of_nodes()
    print('Cantidad inicial de nodos de la red %s:\t%d'%(s,N0))
    
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component=giant.number_of_nodes()
    while largest_component>1:
        D=dict(nx.betweenness_centrality(red))
        maxDegreeNode = max(D, key=D.get) # se toma el nodo de máxima centralidad.
        red.remove_node(maxDegreeNode)
        giant = max(nx.connected_component_subgraphs(red), key=len)
        largest_component=giant.number_of_nodes()
        N=red.number_of_nodes()
        
        print('\rNodos restantes:\t%4d'%(N),end='')
        x_shortest_path.append(1-N/N0)
        y_shortest_path.append(largest_component/N)
    
    '''
    ######################################################
    # Current flow:
    print('\nCURRENT FLOW:')
    red=redes[s].copy()
    x_current_flow=[]
    y_current_flow=[]
    N0=red.number_of_nodes()
    print('Cantidad inicial de nodos de la red %s:\t%d'%(s,N0))
    
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component=giant.number_of_nodes()
    while largest_component>1:
        D=dict(nx.current_flow_closeness_centrality(red))
        maxDegreeNode = max(D, key=D.get) # se toma el nodo de máxima centralidad.
        red.remove_node(maxDegreeNode)
        giant = max(nx.connected_component_subgraphs(red), key=len)
        largest_component=giant.number_of_nodes()
        N=red.number_of_nodes()
        
        print('\rNodos restantes:\t%4d'%(N),end='')
        x_current_flow.append(1-N/N0)
        y_current_flow.append(largest_component/N)
    
    
    '''

#    x_essential=[]
#    y_essential=[]

    ############################################################################
    ############################################################################
    # Gráfico:
    plt.figure()
    plt.plot(x_random,y_random,label='Random')
    plt.plot(x_degree,y_degree,label='Degree')
    plt.plot(x_eigenvector,y_eigenvector,label='Eigenvector')
    plt.plot(x_subgraph,y_subgraph,label='Subgraph')
    plt.plot(x_shortest_path,y_shortest_path,label='Shortest path')
#    plt.plot(x_current_flow,y_current_flow,label='Current flow')
#    plt.plot(x_essential,y_essential,label='Shortest path')
    plt.xticks(np.arange(0, 1.1, step=0.1))
    plt.yticks(np.arange(0, 1.1, step=0.1))
    plt.xlim((0, 0.35))
    plt.tick_params(axis='both', which='major', labelsize=NumberSize)
    plt.xlabel('Fraction of nodes', fontsize=AxisLabelSize)
    plt.ylabel('Largest connected component', fontsize=AxisLabelSize)   
    plt.title('Red %s'%(s), fontsize=TitleSize)
    plt.legend(loc='upper right', fontsize=LegendSize) 
    plt.grid(axis='both', color='k', linestyle='dashed', linewidth=2, alpha=0.1)
    plt.show()
#degree_centrality(G)
#closeness_centrality(G[, u, distance, ...])
#eigenvector_centrality(G[, max_iter, tol, ...])

#%%

#Tabla 3 de Zotenko et. al. (2008)

#dict guardo la fraccion de nodos en la componente gigante luego de eliminar los nodos esenciales de la red 
fraction_nodes = {}

for k in redesStr:
    red_k = redes[k].copy()
    red_k.remove_nodes_from(list(nodoi for nodoi in red_k if nodoi in essentials))
    giant = max(nx.connected_component_subgraphs(red_k), key=len)
    largest_component=giant.number_of_nodes()
    total_nodes = red_k.number_of_nodes()
    fraction_nodes[k]=largest_component/total_nodes

#me guardo las redes sin los nodos esenciales 
redes_ne = {}

#me guardo las redes luego de haberle quitado los nodos rnd no esenciales
redes_rnd = {}

#me guardo los nodos de las redes sin nodos esenciales y sus grados
redes_ne_deg = {}

#me guardo los grados de los nodos esenciales
essentials_deg = {}

#cuento los nodos de las redes esenciales 
degreeCount = {}

#output: fraccion de nodos en la componente gigante luego de sacar los nodos random nonessential
fraction_nodes_rnd = {}

def group_by_degree(red):
    newlist, dicpos = [],{}
    for val, j in red:
        if j in dicpos:
            newlist[dicpos[j]].append(val)
        else:
            newlist.append([val])
            dicpos[j] = len(dicpos)
    return newlist, dicpos


for k in redesStr:
    print(k)
    redes_rnd[k] = redes[k].copy()
    redes_ne[k] = redes[k].copy()
    redes_ne[k].remove_nodes_from(list(nodoi for nodoi in redes_ne[k] if nodoi in essentials))
    redes_ne_deg[k] = sorted(redes_ne[k].degree, key=lambda x: x[1])
    essentials_deg[k] = sorted(list(redes[k].degree(nodoi) for nodoi in essentials if nodoi in redes[k]))
    degreeCount[k] = collections.Counter(essentials_deg[k])
    gnodes, dicpos = group_by_degree(redes_ne_deg[k])
    chosen_nodes = []
    for j in degreeCount[k].keys():
        if dicpos.get(j) is not None:
            print(degreeCount[k][j])
            if len(gnodes[dicpos[j]])>=degreeCount[k][j]:
                chosen_nodes.append(sample(gnodes[dicpos[j]], degreeCount[k][j]))
            else:
                print('alt')
                chosen_nodes.append(sample(gnodes[dicpos[j]], len(gnodes[dicpos[j]])))
    redes_rnd[k].remove_nodes_from([item for sublist in chosen_nodes for item in sublist])
    giant = max(nx.connected_component_subgraphs(redes_rnd[k]), key=len)
    largest_component=giant.number_of_nodes()
    total_nodes = redes_rnd[k].number_of_nodes()
    fraction_nodes_rnd[k]=largest_component/total_nodes
            
