#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 14:58:22 2018

@author: heli
"""
#%%

#librerias

import pandas as pd
import itertools as itr
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import matplotlib.patches as mpatches

#%%

#cargo las tablas con los TCRs asociados a epitopes conocidos de la base de datos publica VDJdb
#filtro por las secuencias con Score mayor a 0 (con algun grado de confianza para su anotacion), Gene igual a TRB == T-cell Receptor Beta, Species igual a Homo Sapiens y MHC_class (clase del Complejo Mayor de Histocompatibilidad) igual a I. 
#me quedo con las columnas CDR3 (Complementary Determining Region 3) que corresponde a la secuencia de aminoacidos del TRB, Epitope que corresponde a la secuencia del epitope y Epitope_gene que corresponde al gen del cual proviene dicho epitope.

path = '/home/heli/Documents/DataEpitopes/'

virusStr = ['YFV', 'CMV', 'EBV', 'HIV-1', 'InfAV', 'HCV']
epitopes = {}

for i in virusStr:    
    df0 = pd.read_csv(path + i, sep='\t')
    df1 = df0[(df0.Score > 0) &(df0.Gene == 'TRB') & (df0.Species == 'HomoSapiens') & (df0.MHC_class == 'MHCI')]
    df2 = df1 [['CDR3', 'Epitope', 'Epitope_gene', 'Epitope_species']]
    df3 = df2.drop_duplicates(subset = 'CDR3')
    epitopes[i] = df3

#%%% 
    
#prueba de pairwise2 de Biopython
#pairwise sequence alignment

alignments = pairwise2.align.localds("LSPADKTNVKAA", "PEEKSAV", blosum62, -10, -1)

#local -> uses Smith Waterman local alignment algorithm
#match parameter 'd': A dictionary returns the score of any pair of characters -> BLOSUM62
#gap parameter 's': Same open and extend gap penalties for both sequences -> in this case, gap open penalty of 10 and a gap extension penalty of 1

#score_only=True prints only the best alignment score

#investigar que matriz BLOSUM usar para alignment de los CDR3 y que parametros para gap penalty

#%%%

#tomo 145 TRB de cada categoria (virus) de forma random 
#genero un dataframe conjunto con 1200 TRB

epitopes_rnd0 = []

for i in virusStr:
    epitopes_rnd0.append(epitopes[i].sample(n=145))

epitopes_rnd = pd.concat(epitopes_rnd0)

epitopes_rnd = epitopes_rnd.drop_duplicates(subset = 'CDR3')
#hay un solo CDR3 que esta repetido, ver!

epitopes_rnd = epitopes_rnd.reset_index(drop=True)

#tomo todos los pares unicos de secuencias y hago pairwise alignment de las mismas para obtener una score que refleje la similud de secuencia

CDR3_pairs = itr.combinations(epitopes_rnd["CDR3"], 2)
seq_scores = []
selected_CDR3_pairs = []

for (seq1, seq2) in CDR3_pairs:
    score = pairwise2.align.localds(seq1, seq2, blosum62, -10, -1, score_only=True)
    seq_scores.append(score)
    if score > 45:
        selected_CDR3_pairs.append((seq1, seq2, {'weight': score}))

max(seq_scores)
min(seq_scores)
np.median(seq_scores)
np.mean(seq_scores)
np.std(seq_scores)

#hay que evaluar adonde poner el cutoff para el score calculado
#por el momento lo fije en 45 que es mas que media + 1 desvio std

#hago un histograma con los scores de las secuencias

plt.figure()
plt.hist(seq_scores, bins=500)
plt.xticks(np.arange(0, 120, step=5))
plt.tick_params(axis='both', which='major', labelsize=10)
plt.xlabel('Sequence score', fontsize=20)
plt.show()

#armo el grafo con el output del pairwise alignment de secuencias

TRB = nx.Graph()
TRB.add_edges_from(selected_CDR3_pairs)

for n in epitopes_rnd['CDR3']:
    if n in list(TRB.nodes()):
        TRB.nodes[n]['epitope_sp'] = epitopes_rnd['Epitope_species'][epitopes_rnd['CDR3'] == n].to_string(index=False)
    else:
        TRB.add_node(n, epitope_sp = epitopes_rnd['Epitope_species'][epitopes_rnd['CDR3'] == n].to_string(index=False))


#len(np.unique(epitopes_rnd['CDR3']))
#len(epitopes_rnd['CDR3'])

TRB.number_of_nodes()
TRB.number_of_edges()

netDeg = np.array(list(TRB.degree()))
CDR3_netDeg = [int(x) for x in netDeg[:,1]]

np.median(CDR3_netDeg)
np.mean(CDR3_netDeg)
np.std(CDR3_netDeg)

plt.figure()
plt.hist(CDR3_netDeg, bins=60)
plt.xticks(np.arange(0, 350, step=10))
plt.tick_params(axis='both', which='major', labelsize=10)
plt.xlabel('Selected sequence degree', fontsize=20)
plt.show()

# use one of the edge properties to control line thickness
edgewidth = [ d['weight']/20 for (u,v,d) in TRB.edges(data=True)]

plt.figure()
nx.draw(TRB, width = edgewidth, node_size=80, font_size=20)
plt.suptitle('TRB sequence similarity network', fontsize=20)
plt.show()

color_code = {'YellowFeverVirus':'m', 'CMV':'g', 'HIV-1':'k',  'EBV':'b', 'InfluenzaA':'r', 'HCV':'c'}

TRB_epitopes = list(nx.get_node_attributes(TRB, "epitope_sp").values())

TRB_node_color = []
for i in TRB_epitopes:
    TRB_node_color.append(color_code[i])


pos = nx.kamada_kawai_layout(TRB)

plt.figure()
nx.draw(TRB,
        width=0.1,
        edge_color = 'k',
        node_color = TRB_node_color, 
        node_size=50,
        font_size=10,
        with_labels=False,
        )
plt.show()

plt.figure()
nx.draw(TRB,
        pos,
        width=0.1,
        edge_color = 'k',
        node_color = TRB_node_color, 
        node_size=50,
        font_size=10,
        with_labels=False,
        )
YFV = mpatches.Patch(color='m', label='Yellow Fever Virus')
CMV = mpatches.Patch(color='g', label='Citomegalovirus')
HIV = mpatches.Patch(color='k', label='HIV-1')
EBV = mpatches.Patch(color='b', label='Epstein-Barr Virus')
InfluenzaA = mpatches.Patch(color='r', label='Influenza A Virus')
HCV = mpatches.Patch(color='c', label='Hepatitis C Virus')
plt.legend(handles=[YFV, CMV, HIV,  EBV, InfluenzaA, HCV], fontsize=20)
plt.show()










