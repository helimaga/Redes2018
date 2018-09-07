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
pathSanti = '/home/santiago/Documentos/RC/tc01_data/'

path = pathSanti

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

#%% Ejercicio 1

'''
1) Considere las tres redes de interacción de proteínas relevadas para levadura disponibles en la
página de la materia. Se trata de: una red de interacciones binarias (yeast_Y2H.txt), de copertenencia
a complejos proteicos (yeast_AP-MS.txt) y obtenida de literatura (yeast_LIT.txt)
obtenidas del Yeast Interactome Database.
'''
#%% Ejercicio 1: Cargar datos

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
    # len es el numero de nodos
	df1b.ix[s,'Diameter'] = nx.diameter(giant)
    
#a = nx.connected_component_subgraphs(redes[s])
#nx.diameter(giant)
#%% Ej. 1(c)

#%% Ej. 1(d)

#%% Ej. 2

dolphins = nx.read_gml(path + 'dolphins.gml')
dolphinsGender = np.array(ldata(path + 'dolphinsGender.txt'))

for n,g in zip(dolphins,range(0,len(dolphinsGender))):
    dolphins.nodes[n]["gender"] = dolphinsGender[g,1]

#%% 2a

nx.draw(dolphins, 
        width=5, 
        node_color=["blue" if g=="m" else ("red" if g=="f" else "yellow") for g in nx.get_node_attributes(dolphins, "gender").values()], 
        node_size=400,
        font_size=20,
        with_labels=True
       )
plt.show()



#%% Ej. 3

'''
1) Considere la red as-22july06.gml creada por Mark Newman que contiene la estructura de los
sistemas autónomos de internet relevada a mediados de 2006.
'''

Newman = nx.read_gml(path + 'as-22july06.gml')

#%% 3a
'''
a) Encuentre gráficamente la distribución de grado P k como función de k explorando
diferentes alternativas: un bineado lineal o logarítmico, utilizando escalas logarítmicas o
lineales en uno o ambos ejes. Discuta que alternativa permite apreciar mejor el carácter
libre de escala de dicha distribución.
'''
degree_sequence = sorted([d for n, d in Newman.degree()], reverse=True)  # Se arma la sucesión de los grados (ordenada).
degreeCount = collections.Counter(degree_sequence) # Se cuenta cuántos nodos hay con cada grado.
deg, cnt = zip(*degreeCount.items()) # Se almacenan los grados y la cantidad de nodos con ese grado.


fig, ax = plt.subplots()
plt.bar(deg, cnt, color='b')
plt.title("Bineado lineal")
plt.ylabel("Cantidad de nodos")
plt.xlabel("Grado")
plt.grid()


fig, ax = plt.subplots()
plt.bar(deg, cnt, color='b')
ax.set_yscale('log')
plt.title("Bineado logari'tmico")
plt.ylabel("Cantidad de nodos")
plt.xlabel("Grado")
plt.grid()

# Discutir.

'''
b. Utilizando funcionalidad de la librería igraph, estime el exponente de dicha distribución.
'''

#%% Ej. 4

'''
4) Asortatividad
'''

'''
a) Considere la red de colaboraciones científicas (netscience.gml) y la red de internet (as-
july06.gml). Analice si nodos de alto grado tienden a conectarse con nodos de alto grado
o por el contrario suelen conectarse a nodos de bajo grado? (i.e la red es asortativa o
disortativa respecto al grado?).
'''


'''
i. Determine, para nodos de grado k, cuánto vale en media el grado de sus vecinos.
[hint R: se puede estimar primero el grado medio de los vecinos de cada nodo de
la red y luego utilizar aggregate sobre esos datos, que permite estimar cantidades
sobre subconjuntos determinados de datos de acuerdo a diferentes criterios]
'''

'''
ii. Analizar la tendencia observada en un gráfico que consigne dicho valor k nn (k)
como función del grado.

'''

'''
iii.Asumiendo que k_nn (k) = ak^μ , estime el exponente de correlación a partir de
realizar una regresión de log k_nn ~ log k. Asegurese de graficar el fiteo en el
grafico anterior. [hint R: lm permite hacer regresiones lineales]
'''

'''
iv. Considere la red de colaboraciones y la de internet nuevamente Encuentre
cuantitativamente la asortatividad de la red utilizando ahora el estimador
propuesto por Newman (ver PDF).
Para ello tenga encuenta lo desarrollado en las eqs [8.26 – 8.29] del libro de
Newman.Como se corresponde este coeficiente con el estimado en el punto
anterior? A qué se debe?
'''

#%%
'''
b) Corra el script de cálculo (puntos i-iii) para las redes Y2H y AP-MS. Puede explicar lo
que observa en cuanto a la asortatividad reportada?
'''





