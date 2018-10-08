#TRABAJO COMPUTACIONAL 2#

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools

'''
import scipy as sp
from scipy import stats
from matplotlib import pylab
from sklearn.linear_model import LinearRegression
import random
import matplotlib.patches as mpatches
from statsmodels.stats.proportion import proportions_ztest
import collections
import igraph

from matplotlib_venn import venn3, venn3_circles
'''

pathHeli = '/home/heli/Documents/Redes/Practicas/TC_02/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejasBiologicos/tc02Data/'
pathSanti = '/home/santiago/Documentos/RC/tc02Data/'
pathDocente = '?'

path = pathHeli

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

giant = max(nx.connected_component_subgraphs(redes[s]), key=len)
largest_component=giant.number_of_nodes()

x_random=[]
y_random=[]

x_degree=[]
y_degree=[]

x_eigenvector=[]
y_eigenvector=[]

x_subgraph=[]
y_subgraph=[]

x_shortest_path=[]
y_shortest_path=[]

x_current_flow=[]
y_current_flow=[]

x_essential=[]
y_essential=[]

N0=redes[s].number_of_nodes()
print(N0)

while largest_component>1:
    randomNode=np.random.choice(list(redes[s].nodes))
    redes[s].remove_node(randomNode)
    giant = max(nx.connected_component_subgraphs(redes[s]), key=len)
    largest_component=giant.number_of_nodes()
    N=redes[s].number_of_nodes()
    print(N)
    x_random.append(1-N/N0)
    y_random.append(largest_component/N)

plt.plot(x_random,y_random)
#degree_centrality(G)
#closeness_centrality(G[, u, distance, ...])
#eigenvector_centrality(G[, max_iter, tol, ...])

#%%

#Tabla 3 de Zotenko et. al. (2008)






















