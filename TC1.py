'''
Trabajo Computacional 1: Ejercicio 1

'''
#%%

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np


plt.close('all')

def ldata(archivo):
    # función de lectura de las tablas
    f=open(archivo)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)	
    return data

path='/home/heli/Documents/Redes/Practicas/TC_01/' # en caso de que los archivos estén en otra carpeta.

# Se cargan los datos en tres listas:
ListaY2H=ldata(path+'yeast_Y2H.txt')
ListaAPMS=ldata(path+'yeast_AP-MS.txt')
ListaLIT=ldata(path+'yeast_LIT.txt')

# Se crean los grafos a partir de las listas:
Y2H = nx.Graph(ListaY2H)
APMS = nx.Graph(ListaAPMS)
LIT = nx.Graph(ListaLIT)

redes = [Y2H,APMS,LIT]
redesStr = ['Y2H','APMS','LIT']
# Ejercicio 1 (a): Presentación de una comparación gráfica de las tres redes.

plt.figure()

for r in redes:
    nx.draw(r, with_labels=True, font_weight='bold')

plt.show()

# Ejercicio 1 (b): Caracteristicas de las redes

# 1.1 numero total de nodos
for r in redes:
    r.number_of_nodes()

# Y2H.nodes()

# 1.2 numero total de enlaces
for r in redes:
    r.number_of_edges()

#%% 1.b.4: PRUEBA grado medio de la red FALTA para las otras 2
'''
degRedes={}
degredesMean={}

for s, r in zip(redesStr,redes):

    degredes{s}=list(r.degree())
    degredes{s} = np.array(degredes{s})
    degredes{s} = degredes[:,1]
    degredes{s} = degredes{s}.astype(int)
    degredesMean{s} = np.mean(degredes{s})
'''
#%% 1b: grado max/min de la red FALTA para las otras 2
for r in redes:
    degY2HMax = np.max(degY2H)
    degY2HMin = np.min(degY2H)

#%% 2 Dirigida o no dirigida
# RTA: la y2h seria no dirigida. 1 no hay pares repetidos. 
# 2 conocimiento previo sobre el experimento
ListaY2HSort = np.sort(ListaY2H)
ListaY2HSortUnique = np.unique(ListaY2HSort,axis=0)

#%% 5 densidad de la red FALTA para los OTROS!!!

#Y2H.is_directed() como sabee?? Tiro false =)

nx.density(Y2H)

# densidad de grafos tipo
nx.density(nx.empty_graph(10))
nx.density(nx.complete_graph(10))
nx.density(nx.star_graph(10))

#%% 6 Clustering: C triangulo (GLOBAL)

# T=3#triangles#triads

nx.transitivity(Y2H)

#%% 6 Clustering: C triangulo (GLOBAL)

# T=3#triangles#triads
f = nx.clustering(Y2H)

nx.average_clustering(Y2H) # corroborar que es equiv a tomar promedio de f

nx.diameter(nx.star_graph(4))




























































