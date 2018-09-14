'''
Trabajo Computacional 1

'''
#%%

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
from matplotlib import pylab
from sklearn.linear_model import LinearRegression

pathHeli = '/home/heli/Documents/Redes/Practicas/TC_01/' # en caso de que los archivos estén en otra carpeta.
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejas/TPsGrupales/tc01_data/'
pathSanti = '/home/santiago/Documentos/RC/tc01_data/'

path = pathJuancho

plt.close('all')

#%%
# Configuraciones para los gráficos:
plt.rc('text', usetex=False)
plt.rc('font', family='serif')

TitleSize=20
AxisLabelSize=20
LegendSize=11
NumberSize=12
LabelSize=8 # etiquetas de los nodos
NodeSize=50 # tamaño de los nodos

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
# Ejercicio 1: Cargar datos

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
	nx.draw(redes[s], with_labels=True, font_weight='bold',font_size=LabelSize,node_size=NodeSize)
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

'''
2) Considere la red social de 62 delfines de Nueva Zelanda
'''

dolphins = nx.read_gml(path + 'dolphins.gml')
dolphinsGender = np.char.lower(np.array(ldata(path+'dolphinsGender.txt')))

dolphinsSIndex = dolphinsGender[:, 0].argsort()
dolphinsGenderSort=dolphinsGender[dolphinsSIndex]

for n,g in zip(dolphins,range(len(dolphinsGenderSort))):
    dolphins.nodes[n]["gender"] = dolphinsGenderSort[g, 1]
 
dolphins.nodes("gender")


#%% 2a

'''
Examine las diferentes opciones de layout para este grafo e identifique la que le resulte
mas informativa. Justifique su eleccion detallando las caracteristicas estructurales
de la red que su eleccion pone en evidencia. Incluya en la representacion grafica de la
red informacion sobre el sexo de los delfines. 
'''
plt.figure()
nx.draw(dolphins, 
        width=5, 
        node_color=["blue" if g=="m" else ("red" if g=="f" else "yellow") for g in nx.get_node_attributes(dolphins, "gender").values()], 
        node_size=400,
        font_size=20,
        with_labels=True
       )
plt.show()

#%% 2b

'''
Se trata de una red en la cual prevalece la homofilia en la variable genero?
Para responder
i.Considere la distribucion nula para la fraccion de enlaces que vinculan generos diferentes, generada a partir de al menos 1000 asignaciones aleatorias de genero
ii.A partir de lo obtenido proponga una estimacion para el valor y el error de dichas cantidad cuando no existe vinculo entre topologia de la red y asignacion de genero
Compare su estimacion con el valor medio esperado.
iii.Estime la significancia estadistica del valor observado en el caso de la red real
'''
#acá hay un temita con respecto a cómo considerar a los 'na', por el momento si hay un nodo 'f' o'm' que se conecta con un 'na' considero que es un enlace de generos diferentes
#primero cuento los enlaces entre géneros diferentes

edgesxyGender=0

for (u,v) in dolphins.edges():
    if dolphins.node[u]['gender'] != dolphins.node[v]['gender']:
        print (u,v)
        edgesxyGender+=1

nodesTot=dolphins.number_of_nodes()
edgesTot=dolphins.number_of_edges()

#aca calculo la fraccion de enlaces entre generos diferentes para la red de estudio
ratioEdgesGender=edgesxyGender/edgesTot

#ESTRATEGIA 1: PERMUTACION 

#luego hago 1000 asignaciones al azar del género a la red dolphinsRnd
dolphinsRnd=dolphins.copy()
ratioEdgesRndGender=[]
number_of_iter=10000

for i in range(number_of_iter):
    np.random.shuffle(dolphinsGenderSort[:,1])
    for n,g in zip(dolphinsRnd,range(len(dolphinsGenderSort))):
        dolphinsRnd.nodes[n]["gender"] = dolphinsGenderSort[g, 1]
    edgesxyRndGender=0
    for (u,v) in dolphinsRnd.edges():
        if dolphinsRnd.node[u]['gender'] != dolphinsRnd.node[v]['gender']:
            edgesxyRndGender+=1
    ratioEdgesRndGender.append(edgesxyRndGender/edgesTot)

#ploteo la distribucion nula (random) para la fraccion de enlaces entre generos diferentes
plt.figure()
plt.hist(ratioEdgesRndGender, bins=40)
plt.title('',fontsize=TitleSize)
plt.ylabel('Frecuencia',fontsize=AxisLabelSize)
plt.xlabel('Proporcion de enlaces que vinculan generos diferentes',fontsize=AxisLabelSize)
#plt.grid()
plt.tight_layout()
plt.show()

#calculo la probabilidad de que la proporcion de enlaces entre generos distintos sea tan o mas extrema que la observada para la red de estudio considerando que HO es verdadera (considerando como distribucion nula a la generada por permutar los sexos de los delfines, dejando inalterada la topologia de la red)

p_value=sum(i >= ratioEdgesGender for i in ratioEdgesRndGender)/number_of_iter
p_value

#este p-valor estaria indicando que la red tiene una proporcion de enlaces entre generos distintos mucho menor que lo esperado por azar 
#por ende, es homofilica!

#calculo la media y el desvio estadar para la fraccion de enlaces de la dist nula
pop_mu=np.mean(ratioEdgesRndGender)
pop_std=np.std(ratioEdgesRndGender)

#calculo la fraccion de nodos 'f', 'm' y 'na' para la red de estudio
probF=(dolphinsGender[:,1]=='f').sum()/nodesTot
probM=(dolphinsGender[:,1]=='m').sum()/nodesTot
pronNA=(dolphinsGender[:,1]=='na').sum()/nodesTot

#esto siguiente valdria si no tuviesemos la categoria 'na'
#estimo la media y la var de la fraccion de enlaces entre generos diferentes en una red en la cual no tengo en cuenta la topologia de la misma (propuesto en las diapos clase)
true_mu=2*probF*probM
true_var=true_mu*(1-true_mu)

#hago un test-z para proporciones para obtener una significancia estadisticas
from statsmodels.stats.proportion import proportions_ztest

proportions_ztest(edgesxyGender, edgesTot, true_mu, 'two-sided', true_var)
   
sm.qqplot(np.asarray(ratioEdgesRndGender), line='45')
pylab.show()
sp.stats.shapiro(np.asarray(ratioEdgesRndGender))

#muy lejos de ser normal la distribucion de la fraccion de enlaces que unen generos distintos
#por ende, es cualquiera hacer un test-t

# =============================================================================

#G = nx.Graph()
#G.add_node(1,color='red')
#G.add_node(2,color='red')
#G.add_node(3,color='blue')
#G.add_node(4,color='blue')
#G.add_node(5,color='red')
#
#
#G.add_edges_from([(1,2),(1,3),(3,4),(4,5)])
#
#for (u,v) in G.edges:
#    if G.node[u]['color'] != G.node[v]['color']:
#        print (u,v)
# =============================================================================

'''
#los datos de nuestra muestra/red de estudio, los '1's representan enlaces entre generos diferentes y los '0's entre el mismo genero
datadolphins=np.asarray([1]*63+[0]*96)

#test t de student, prueba a 2 colas 
#haciendo aqui el test-t de student asumo que la media de la fraccion de enlaces entre generos diferentes sigue una distribucion normal, lo cual podria no ser cierto

ttest1 = sp.stats.ttest_1samp(datadolphins, pop_mu)
#rechazaria HO
ttest2 = sp.stats.ttest_1samp(datadolphins, true_mu)
#no rechazaria HO
'''


#%%

# Identifique alguna metodología basada en observables topológicos para eliminar
# nodos secuencialmente de la red de manera de dividirla en dos componentes de tamaños
# comparables en el menor número de pasos. Explique y muestre los resultados obtenidos.
# Intente cuantificar su estrategia comparándola con lo que se obtendría al eliminar nodos
# de manera aleatoria.

dolphinsSplit = dolphins.copy()
edges = list(dolphins.edges())
for ed in edges:
    edGender = (dolphins.nodes[ed[0]]['gender'],dolphins.nodes[ed[1]]['gender'])
    edGender = sorted(edGender, key=str.lower)
    if (edGender == ['f','m']) | (edGender == ['m','na']):
        dolphinsSplit.remove_edge(*ed)
            
plt.figure()
nx.draw(dolphinsSplit, 
        width=2, 
        node_color=["blue" if g=="m" else ("red" if g=="f" else "yellow") for g in nx.get_node_attributes(dolphins, "gender").values()], 
        node_size=40,
        font_size=5,
        with_labels=True
       )



'''
1) Considere la red as-22july06.gml creada por Mark Newman que contiene la estructura de los sistemas autónomos de internet relevada a mediados de 2006.
'''

Newman = nx.read_gml(path + 'as-22july06.gml')

#%% 3a
'''
a) Encuentre gráficamente la distribución de grado P_k como función de k explorando
diferentes alternativas: un bineado lineal o logarítmico, utilizando escalas logarítmicas o
lineales en uno o ambos ejes. Discuta que alternativa permite apreciar mejor el carácter
libre de escala de dicha distribución.
'''
import collections

degree_sequence = sorted([d for n, d in Newman.degree()], reverse=True)  # Se arma la sucesión de los grados (ordenada).
degreeCount = collections.Counter(degree_sequence) # Se cuenta cuántos nodos hay con cada grado.
deg, cnt = zip(*degreeCount.items()) # Se almacenan los grados y la cantidad de nodos con ese grado.

binsLin = np.linspace(min(deg),max(deg),151)
binsLog = np.logspace(np.log10(min(deg)),np.log10(max(deg)),151)
sp=0
#fig, axes = plt.subplots(2, 2, subplot_kw=dict(polar=True))
#for binning in [binsLog,binsLin]:
#    for scale in ['log','lin']:
#        sp=sp+1
#        fig, ax = plt.subplots(2,2,sp)
#        plt.hist(deg, bins=binning, edgecolor='k')
#        ax.set_xscale(scale)
#        ax.set_yscale(scale)

fig, ax = plt.subplots()
plt.bar(deg, cnt, color='b', edgecolor='k')
ax.set_xscale('log')
ax.set_yscale('log')
plt.tick_params(axis='both', which='major', labelsize=NumberSize)
plt.legend(loc='best',fontsize=LegendSize)
plt.title(r'Bineado lineal',fontsize=TitleSize)
plt.ylabel(r'Cantidad de nodos',fontsize=AxisLabelSize)
plt.xlabel(r'Grado',fontsize=AxisLabelSize)
plt.grid()
plt.tight_layout()
plt.show()

fig, ax = plt.subplots()
plt.bar(deg, cnt, color='b')
ax.set_xscale('log')
ax.set_yscale('log')
plt.tick_params(axis='both', which='major', labelsize=NumberSize)
#plt.legend(loc='best',fontsize=LegendSize)
plt.title(r'Bineado logar\'itmico',fontsize=TitleSize)
plt.ylabel(r'Cantidad de nodos',fontsize=AxisLabelSize)
plt.xlabel(r'Grado',fontsize=AxisLabelSize)
#plt.grid()
plt.tight_layout()
plt.show()

# Discutir.

'''
b. Utilizando funcionalidad de la librería igraph, estime el exponente de dicha distribución.
'''
import igraph

# fit_power_law fits a power-law distribution to a data set. 
#fit_power_law(x, xmin = NULL, start = 2, force.continuous = FALSE, implementation = c("plfit", "R.mle"))

fit = igraph.fit_power_law(deg+1, 1)

#g <- barabasi.game(1000) # increase this number to have a better estimate
#d <- degree(g, mode="in")
#fit1 <- fit_power_law(d+1, 10)

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
#%% Leer redes

redesStr = ['netscience','as-22july06']
redes = {}
avnd = {}
degree = {}
nodes = {}
degreeAvnd = {}

for s in redesStr: # cambiar!!!
	redes[s] = nx.read_gml(path + s + '.gml')

nodes = redes[s].nodes()
avnd = nx.average_neighbor_degree(redes[s])
degree0 = redes[s].degree()

ds = [dict(degree0), avnd]
d = {}
for k in avnd.keys():
    d[k] = list(d[k] for d in ds)

degreeAvnd = d.values()
degreeAvnd = np.array(list(degreeAvnd))

#%% Grafico

# x from 0 to 30

x = np.transpose(degreeAvnd[:,0])
y = np.transpose(degreeAvnd[:,1])

x = x.reshape((-1,1))  # conversion between (N,) dimension arrays and (N,1) 
y = y.reshape((-1,1))

logx = np.log10(x)
logy = np.log10(y)

logx = logx.reshape((-1,1))
logy = logy.reshape((-1,1))

model = LinearRegression()
model.fit(logx, logy)

logx_new = np.linspace(min(logx), max(logx), 100)
logy_new = model.predict(logx_new[:, np.newaxis])

plt.figure()
plt.plot(x, y,'.k')
plt.plot(10**logx_new, 10**logy_new,'r')
plt.xlabel('Degree, k')
plt.ylabel('Avergage degree of neighbors, knn')
plt.xscale('log')
plt.yscale('log')


plt.show()

#%%
