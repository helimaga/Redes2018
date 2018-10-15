###

#TRABAJO COMPUTACIONAL 2#

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
import collections
from random import sample
#import scipy as sp
from sklearn.linear_model import LinearRegression

'''
from scipy import stats
from matplotlib import pylab
import matplotlib.patches as mpatches
from statsmodels.stats.proportion import proportions_ztest
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

#Punto B
#Caracteristicas de las redes analizadas - Tabla 1 de Zotenko et. al. (2008)

dfB1 = pd.DataFrame()

for s in redesStr:
	dfB1.loc[s,'Nodes'] = redes[s].number_of_nodes()
	dfB1.loc[s,'Edges'] = redes[s].number_of_edges()
	if len(redes[s]) == len(np.unique(redes[s],axis=0)):
		dfB1.loc[s,'Directionality'] = 'Prob-Undir'
	else:
		dfB1.loc[s,'Directionality'] = 'Dir'
	netDeg = np.array(list(redes[s].degree()))
	netDeg = netDeg[:,1]
	netDeg = netDeg.astype(int)
	dfB1.loc[s,'DegMean'] = np.mean(netDeg)
	dfB1.loc[s,'DegMin'] = np.min(netDeg)
	dfB1.loc[s,'DegMax'] = np.max(netDeg)
	dfB1.loc[s,'DegDensity'] = nx.density(redes[s])
	dfB1.loc[s,'C_tri'] = nx.transitivity(redes[s])
	dfB1.loc[s,'C_avg'] = nx.average_clustering(redes[s])
	giant = max(nx.connected_component_subgraphs(redes[s]), key=len)
	dfB1.loc[s,'Diameter (Giant Subgraph)'] = nx.diameter(giant)

tabla1 = dfB1.to_latex(buf=None, columns=['Nodes','Edges','DegMean','C_avg'], col_space=None, bold_rows=False,float_format='%.3f')


#%%

#Punto B
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

<<<<<<< HEAD
tabla=dfB2.to_latex(buf=None, columns=['Y2H','AP-MS','LIT','LIT_Reguly'], col_space=None, bold_rows=False,float_format='%.3f')

print('Tabla 2:\n')
print(tabla)
=======
tabla2 = dfB2.to_latex(buf=None, columns=['Y2H','AP-MS','LIT','LIT_Reguly'], col_space=None, bold_rows=False,float_format='%.3f')

>>>>>>> Tabla 3 reformulada

#%% 

#Punto B
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
#plt.title('Figure 1A. Relationship between degree and essentiality in the tested networks', fontsize=15)
plt.legend(loc='upper right', fontsize=15) 
plt.grid(axis='both', color='k', linestyle='dashed', linewidth=2, alpha=0.1)
plt.show()
plt.savefig(path+'/Figuras/Figura1.pdf')


#%%

#Punto C - Analisis de vulnerabilidad 
#Figura 3 de Zotenko et. al. (2008)

def RemoveNodes(RED,method):
    # Función que remueve los nodos.
    red=RED.copy()

    print('\n',method)
    x=[]
    y=[]
    N0=red.number_of_nodes()
    print('Cantidad inicial de nodos de la red %s:\t%d'%(s,N0))
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component=giant.number_of_nodes()
    
    # Se remueven los nodos según el método:
    while largest_component>1:
        if method=='Random':
            Node=np.random.choice(list(red.nodes))
        elif method=='Degree':
            D=dict(nx.degree_centrality(red))
            Node = max(D, key=D.get)
        elif method=='Eigenvector':
            D=dict(nx.eigenvector_centrality(red,tol=1e-03))
            Node = max(D, key=D.get)
        elif method=='Subgraph':
            D=dict(nx.subgraph_centrality(red))
            Node = max(D, key=D.get)
        elif method=='Shortest path':
            D=dict(nx.betweenness_centrality(red))
            Node = max(D, key=D.get)
        elif method=='Current flow':
            D=dict(nx.current_flow_betweenness_centrality(red))
            Node = max(D, key=D.get)
            
        red.remove_node(Node)
        giant = max(nx.connected_component_subgraphs(red), key=len)
        largest_component=giant.number_of_nodes()
        N=red.number_of_nodes()
        print('\rNodos restantes:\t%4d'%(N),end='')
        x.append(1-N/N0)
        y.append(largest_component/N)
    
    np.save(path + '%s-x_%s'%(s,method), x)
    np.save(path + '%s-y_%s'%(s,method), y)
    
    plt.plot(x,y,label=method)
    return

def RemoveEssentials(RED):
    # Función que remueve los nodos esenciales de una vez.
    red=RED.copy()

    N0=red.number_of_nodes()
    print('N0 = ',N0)
    red.remove_nodes_from(essentials)
    N=red.number_of_nodes()
    print('N = ',N)
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component=giant.number_of_nodes()
    print('largest_component = ',largest_component)
    

    x = 1-N/N0
    y = largest_component/N0
    print('x,y = ',x,y)
    plt.plot(x,y,'o',label='Essentials')
    return

for s in redesStr:
    print('\n\nRed: %s'%(s))
    plt.figure()
    
    RemoveEssentials(redes[s])
    
    methods=['Random', 'Degree', 'Eigenvector', 'Subgraph', 'Shortest path']#, 'Current flow']
    
    for m in methods:
        RemoveNodes(redes[s],m)
    
    
    plt.xticks(np.arange(0, 1.1, step=0.1))
    plt.yticks(np.arange(0, 1.1, step=0.1))
    plt.xlim((0, 0.42))
    plt.tick_params(axis='both', which='major', labelsize=NumberSize)
    plt.xlabel('Fraction of nodes', fontsize=AxisLabelSize)
    plt.ylabel('Largest connected component', fontsize=AxisLabelSize)   
    plt.title('Red %s'%(s), fontsize=TitleSize)
    plt.legend(loc='upper right', fontsize=LegendSize)
    plt.grid(axis='both', color='k', linestyle='dashed', linewidth=2, alpha=0.1)
    plt.show()
    plt.savefig(path+'/Figuras/Figura3-%s.pdf'%(s))



#%%

#Punto C - Analisis de vulnerabilidad 
#Tabla 3 de Zotenko et. al. (2008)

def fraction_nodes_giantcomp(red0, chosen_nodes):
    
    red = red0.copy()
    red.remove_nodes_from(chosen_nodes)
    giant = max(nx.connected_component_subgraphs(red), key=len)
    largest_component = giant.number_of_nodes()
    total_nodes = red.number_of_nodes()
    frac_nodes = largest_component/total_nodes
    
    return frac_nodes

#me guardo las redes sin los nodos esenciales 
redes_ne = {}

#me guardo las redes luego de haberle quitado los nodos rnd no esenciales
redes_rnd = {}

#me guardo los nodos de las redes sin nodos esenciales y sus grados
redes_ne_deg = {}

#me guardo los nodos esenciales y sus grados
essentials_deg = {}


def group_by_degree_deciles(node_deg, s):
    nodes_degdec = {}
    j = 10
    m = []
    for n in range(len(node_deg)):
        if node_deg[n][1] <= degree_cutoff[s][j]:
            m.append(node_deg[n][0])
            if n == len(node_deg)-1:
                nodes_degdec[j] = m
        else:
            if degree_cutoff[s][j] == degree_cutoff[s][j+10]:
                p = 0
                while degree_cutoff[s][j+p] == degree_cutoff[s][j+p+10]:
                    p += 10
                j += p 
                nodes_degdec[j] = m
                m = []
                m.append(node_deg[n][0])
                j += 10
            else:
                nodes_degdec[j] = m
                m = []
                m.append(node_deg[n][0])
                j += 10
    return nodes_degdec

dftabla3 = pd.DataFrame()
rnd_iter = 30

for s in redesStr:
    
    redes_ne[s] = redes[s].copy()
    redes_ne[s].remove_nodes_from(list(nodei for nodei in redes_ne[s] if nodei in essentials))
    redes_ne_deg[s] = sorted(redes_ne[s].degree, key=lambda x: x[1])
    essentials_deg[s] = sorted(list((nodei, redes[s].degree(nodei)) for nodei in redes[s] if nodei in essentials), key=lambda x: x[1])
    
    dec_nonEssentials = group_by_degree_deciles(redes_ne_deg[s], s)
    dec_Essentials = group_by_degree_deciles(essentials_deg[s], s)
    
    fraction_rnd_nodes = []
    
    for i in range(rnd_iter):
        
        chosen_rnd_nodes = []
        chosen_nodes = []
        for j in dec_Essentials.keys():
            numEssentials = len(dec_Essentials[j])
            numNonEssentials = len(dec_nonEssentials[j])
            if numEssentials <= numNonEssentials:
                chosen_rnd_nodes.append(sample(dec_nonEssentials[j], numEssentials))
                if i == rnd_iter-1:
                    chosen_nodes.append(dec_Essentials[j])
            else:
                chosen_rnd_nodes.append(sample(dec_nonEssentials[j], numNonEssentials))
                if i == rnd_iter-1:
                    chosen_nodes.append(sample(dec_Essentials[j], numNonEssentials))
              
        chosen_rnd_nodes = [item for sublist in chosen_rnd_nodes for item in sublist]  
        fraction_rnd_nodes.append(fraction_nodes_giantcomp(redes[s], chosen_rnd_nodes))
    
    chosen_nodes = [item for sublist in chosen_nodes for item in sublist] 
     
    dftabla3.loc[s, 'Essential'] = fraction_nodes_giantcomp(redes[s], chosen_nodes)
    dftabla3.loc[s, 'Random nonessential Mean'] = np.mean(fraction_rnd_nodes)
    dftabla3.loc[s, 'Random nonessential Std'] = np.std(fraction_rnd_nodes)
 
tabla3=dftabla3.to_latex(float_format='%.3f')

#%%

#Punto D - Esencialidad
#Figura 2B de He et al. (2006)

redes_deg = {}
essentials_deg = {}
degCount = {}
degCount_essentials = {}
freq_pe = {}

x = {}
y = {}
model = {}
m = {}
b = {}

alpha = {}
beta = {}

for s in redesStr:
    essentials_deg[s] = sorted(list(redes[s].degree(nodoi) for nodoi in redes[s] if nodoi in essentials))
    redes_deg[s] = sorted(list(redes[s].degree(nodoi) for nodoi in redes[s]))
    degCount[s] = collections.Counter(redes_deg[s])
    degCount_essentials[s] = collections.Counter(essentials_deg[s])
    freq_pe[s] = dict()
    for l in degCount[s].keys():
        if degCount_essentials[s].get(l) is not None:
            freq_pe[s][l] = np.log(1 - (degCount_essentials[s][l]/degCount[s][l]))
        else:
            freq_pe[s][l] = 0

    x = np.array(list(freq_pe[s].keys())[:9])
    y = np.array(list(freq_pe[s].values())[:9])
    
    x = x.reshape((-1,1))
    y = y.reshape((-1,1))

    model[s] = LinearRegression()
    model[s].fit(x, y)
    
    y_predict = model[s].predict(x)
    r2 = round(float(model[s].score(x,y)),4)
    
    m[s] = round(float(model[s].coef_),4)
    b[s] = round(float(model[s].intercept_),4)

    alpha[s] = 1 - np.exp(m[s])
    beta[s] =  1 - np.exp(b[s])
    
    textStr= '$ln(1-P_{E})=%.2fk%.2f$\n$alpha=%.2f$\n$beta=%.2f$\n$r^{2}=%.2f$'%(m[s],b[s],alpha[s],beta[s],r2)
    
    plt.figure()
    plt.plot(x, y,'.k')
    plt.plot(x, y_predict, 'r', label=textStr)
    plt.xlabel(r'Degree or protein connectivity ($k$)', fontsize=20)
    plt.ylabel(r'$ln(1-P_{E})$', fontsize=20)
    plt.title('Red ' + s ,fontsize=30)
    plt.legend(fontsize=20)
    plt.show()

#%%

#Punto D - Esencialidad
#Tabla 5 de Zotenko et. al. (2008) - En Número esperado solo deben incluir el obtenido
#a partir del ajuste lineal.

def p_same_pair(alpha, beta, ki, kj):
    p_essential_i = 1-((1-alpha)**ki)*(1-beta)
    p_essential_j = 1-((1-alpha)**kj)*(1-beta)
    return 1 - p_essential_i - p_essential_j + 2*p_essential_i*p_essential_j

total_pairs = {}
same_pairs = {}
expected_same_pairs = {}

for s in redesStr:   
    proba = []
    
    if s == 'Y2H':
        common_neighbors = 1
    else:
        common_neighbors = 3
        
    node_pairs = itertools.combinations(redes[s].nodes, 2)
    
    number_of_pairs = 0
    number_of_epairs = 0
    number_of_nonepairs = 0
    
    for (nodei, nodej) in node_pairs:
        if nodej not in redes[s].neighbors(nodei):
            path_len2 = list(nx.all_simple_paths(redes[s], nodei, nodej, cutoff=2))
            if len(path_len2) >= common_neighbors:
                proba.append(p_same_pair(alpha[s], beta[s], redes[s].degree(nodei), redes[s].degree(nodej)))
                number_of_pairs += 1
                if nodei in essentials:
                    if nodej in essentials: 
                        number_of_epairs +=1
                elif nodej not in essentials:
                    number_of_nonepairs += 1
    
    number_of_same_pairs = number_of_epairs + number_of_nonepairs
    total_pairs[s] = number_of_pairs
    same_pairs[s] = number_of_same_pairs
    expected_same_pairs[s] = np.sum(proba)

#%%

'''
for s in redesStr:
    redes_rnd[s] = redes[s].copy()
    redes_ne[s] = redes[s].copy()
    redes_ne[s].remove_nodes_from(list(nodoi for nodoi in redes_ne[s] if nodoi in essentials))
    redes_ne_deg[s] = sorted(redes_ne[s].degree, key=lambda x: x[1])
    essentials_deg[s] = sorted(list(tuple(nodoi, redes[s].degree(nodoi)) for nodoi in redes[s] if nodoi in essentials))
    degreeCount[s] = collections.Counter(essentials_deg[s])
    gnodes, dicpos = group_by_degree(redes_ne_deg[s])
    chosen_nodes = []
    conflict_keys = []
    for j in degreeCount[s].keys():
        
        numEssentials = degreeCount[s][j]
        
        if dicpos.get(j) is not None: 
            
            numNonEssentials = len(gnodes[dicpos[j]])
            
            if numNonEssentials >= numEssentials:
                chosen_nodes.append(sample(gnodes[dicpos[j]], numEssentials))
            else:
                chosen_nodes.append(sample(gnodes[dicpos[j]], numNonEssentials))
                diff = numEssentials - numNonEssentials
                conflict_keys.append(j)
                conflict_keys.append(j)
        else:
            conflict_keys.append(j)
    print(conflict_keys)
                
    redes_rnd[s].remove_nodes_from([item for sublist in chosen_nodes for item in sublist])
    giant = max(nx.connected_component_subgraphs(redes_rnd[s]), key=len)
    largest_component=giant.number_of_nodes()
    total_nodes = redes_rnd[s].number_of_nodes()
    fraction_nodes_rnd[s]=largest_component/total_nodes

'''


