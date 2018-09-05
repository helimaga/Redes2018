'''
Trabajo Computacional 1: Ejercicio 1

'''
#%%

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

pathHeli = '/home/heli/Documents/Redes/Practicas/TC_01/' # en caso de que los archivos estén en otra carpeta.
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejas/TC01/tc01_data/'
pathSanti = ''

path = pathHeli

plt.close('all')

#%% Ejercicio 1

'''
1) Considere las tres redes de interacción de proteínas relevadas para levadura disponibles en la
página de la materia. Se trata de: una red de interacciones binarias (yeast_Y2H.txt), de copertenencia
a complejos proteicos (yeast_AP-MS.txt) y obtenida de literatura (yeast_LIT.txt)
obtenidas del Yeast Interactome Database.
'''
#%% Ejercicio 1: Cargar datos
def ldata(archivo):
    # función de lectura de las tablas
    f=open(archivo)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)	
    return data

redesStr = ['Y2H','AP-MS','LIT']
redes = {}

for s in redesStr:
	redes[s] = nx.Graph(ldata(path + 'yeast_' + s + '.txt'))


#%% Ej. 1(a)
'''
Presente una comparación gráfica de las 3 redes.
'''

for s in redesStr:
	plt.figure()
	nx.draw(redes[s], with_labels=True, font_weight='bold')
	plt.show()
#%% Ej. 1(b)
'''
1b. Resuma en una tabla las siguientes características de dichas redes
i. El número total de nodos, N
ii. El número total de enlaces L, de la red
iii. Si se trata de una red dirigida o no-dirigida
iv. El grado medio/máximo/mínimo de la red
v. La densidad de la red
vi. Los coeficientes de clustering <Ci> y C_delta de la red.
vii. Diámetro de la red.
'''

df1b = pd.DataFrame()

for s in redesStr:
	df1b.ix[s,'Nodes'] = redes[s].number_of_nodes()
	df1b.ix[s,'Edges'] = redes[s].number_of_edges()
	if len(redes[s]) == len(np.unique(redes[s],axis=0)):
		df1b.ix[s,'Directionality'] = 'Prob-Undir'
	else:
		df1b.ix[s,'Directionality'] = 'Dir'
	netDeg = np.array(list(redes[s].degree()))
	netDeg = netDeg[:,1]
	netDeg = netDeg.astype(int)
	df1b.ix[s,'DegMean'] = np.mean(netDeg)
	df1b.ix[s,'DegMin'] = np.min(netDeg)
	df1b.ix[s,'DegMax'] = np.max(netDeg)
	df1b.ix[s,'DegDensity'] = nx.density(redes[s])
	df1b.ix[s,'ClustGlob'] = nx.transitivity(redes[s])
	df1b.ix[s,'ClustLoc'] = nx.average_clustering(redes[s]) # corroborar que es equiv a tomar promedio de f
	giant = max(nx.connected_component_subgraphs(redes[s]), key=len)	
	df1b.ix[s,'Diameter'] = nx.diameter(giant)
    
    