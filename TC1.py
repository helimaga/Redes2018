'''
Trabajo Computacional TC_01

Heli Magali Garcia Alvarez
Santiago Scheiner
Juan Ignacio Gossn
'''
#%% Bloque de inicializacion
# Importamos todos los modulos que necesitaremos para el trabajo:

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
from matplotlib import pylab
from sklearn.linear_model import LinearRegression
import random
import matplotlib.patches as mpatches
from statsmodels.stats.proportion import proportions_ztest
import collections
import igraph
import itertools

# Seleccion de path según la máquina de trabajo:

pathHeli = '/home/heli/Documents/Redes/Practicas/TC_01/'
pathJuancho = '/home/gossn/Dropbox/Documents/Materias_doctorado/RedesComplejas/TPsGrupales/tc01_data/'
pathSanti = '/home/santiago/Documentos/RC/tc01_data/'
pathDocente = '?'

path = pathSanti



# Configuraciones para los graficos:
plt.close('all')
plt.rc('text', usetex=False)
plt.rc('font', family='serif')

TitleSize=20
AxisLabelSize=20
LegendSize=11
NumberSize=12
LabelSize=8 # etiquetas de los nodos
NodeSize=50 # tamaño de los nodos

# Definición de función de lectura de las tablas:

def ldata(archivo):

    f=open(archivo)
    data=[]
    for line in f:
        line=line.strip()
        col=line.split()
        data.append(col)	
    return data

#%% Ej. 1

'''
1) Considere las tres redes de interacción de proteínas relevadas para levadura disponibles en la
página de la materia. Se trata de: una red de interacciones binarias (yeast_Y2H.txt), de copertenencia
a complejos proteicos (yeast_AP-MS.txt) y obtenida de literatura (yeast_LIT.txt)
obtenidas del Yeast Interactome Database.
'''
# Cargar datos

redesStr = ['Y2H','AP-MS','LIT']
redes = {}

for s in redesStr:
	redes[s] = nx.Graph(ldata(path + 'yeast_' + s + '.txt'))


#%% Ej. 1(a)
'''
Presente una comparación gráfica de las 3 redes.
'''

fig, axs = plt.subplots(1,len(redesStr))
axs = axs.ravel()
sp = -1
for s in redesStr:
    sp += 1

    plt.sca(axs[sp])
    # draw all nodes homogeneously, and edge weights as filtered
    nx.draw(redes[s], with_labels=False, node_size=3, axs=axs[sp])

    axs[sp].set_title(s, fontsize=10)
    axs[sp].set_axis_off()

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

# Armamos un dataframe, correspondiente al modulo "pandas"
df1b = pd.DataFrame()

for s in redesStr:
    # i
	df1b.loc[s,'Nodes'] = redes[s].number_of_nodes()
    # ii
	df1b.loc[s,'Edges'] = redes[s].number_of_edges()
	# iii
	if len(redes[s]) == len(np.unique(redes[s],axis=0)):
		df1b.loc[s,'Directionality'] = 'Prob-Undir'
	else:
		df1b.loc[s,'Directionality'] = 'Dir'
    # iv
	netDeg = np.array(list(redes[s].degree()))
	netDeg = netDeg[:,1]
	netDeg = netDeg.astype(int)
	df1b.loc[s,'DegMean'] = np.mean(netDeg)
	df1b.loc[s,'DegMin'] = np.min(netDeg)
	df1b.loc[s,'DegMax'] = np.max(netDeg)
    # v
	df1b.loc[s,'DegDensity'] = nx.density(redes[s])
    # vi
	df1b.loc[s,'ClustGlob'] = nx.transitivity(redes[s])
	df1b.loc[s,'ClustLoc'] = nx.average_clustering(redes[s])
    # vii Dado que las redes son disconexas, estimamos el diametro de la componente gigante
	giant = max(nx.connected_component_subgraphs(redes[s]), key=len)
	df1b.loc[s,'Diameter'] = nx.diameter(giant)

print(df1b)

#%% Ej. 1(c)
'''
Teniendo en cuenta la naturaleza de las interacciones reportadas, diga si es razonable lo
que encuentra para ciertos observables calculados.
'''
#%% Ej. 1(d)
'''
Construya un diagrama de Venn que permita reconocer la cobertura, especificidad y
coherencia de las interacciones reportadas por los tres datasets
'''
#%% Ej. 2
'''
Considere la red social de 62 delfines de Nueva Zelanda
'''
# Cargar datos

dolphins = nx.read_gml(path + 'dolphins.gml') # red
dolphinsGender = np.array(ldata(path+'dolphinsGender.txt')) # numpy array con los generos

# En este loop, identificamos el genero de cada delfin nodo por nodo:
for n in dolphins.nodes():
    dolphins.nodes[n]["gender"] = dolphinsGender[dolphinsGender[:,0]==n,1][0]

# Para corroborar, hacemos un print del cada nodo con su respectivo genero
print(dolphins.nodes("gender"))

#numero de nodos y enlaces de la red dolphins
nodesTot=dolphins.number_of_nodes()
edgesTot=dolphins.number_of_edges()


#%% Ej. 2(a)
'''
Examine las diferentes opciones de layout para este grafo e identifique la que le resulte
mas informativa. Justifique su eleccion detallando las caracteristicas estructurales
de la red que su eleccion pone en evidencia. Incluya en la representacion grafica de la
red informacion sobre el sexo de los delfines. 
'''

plt.figure()
nx.draw(dolphins,
        width=1,
        edge_color = 'c',
        node_color=["blue" if g=="m" else ("red" if g=="f" else "yellow") for g in nx.get_node_attributes(dolphins, "gender").values()], 
        node_size=80,
        font_size=20,
        with_labels=True,
       )
male = mpatches.Patch(color='b', label='Male')
female = mpatches.Patch(color='r', label='Female')
unknown = mpatches.Patch(color='y', label='Unknown')
plt.legend(handles=[male,female,unknown])
plt.suptitle('Red Delfines')
plt.show()

#%% Ej. 2(b)

'''
Se trata de una red en la cual prevalece la homofilia en la variable genero?
Para responder
i.Considere la distribucion nula para la fraccion de enlaces que vinculan generos diferentes, generada a partir de al menos 1000 asignaciones aleatorias de genero
ii.A partir de lo obtenido proponga una estimacion para el valor y el error de dichas cantidad cuando no existe vinculo entre topologia de la red y asignacion de genero
Compare su estimacion con el valor medio esperado.
iii.Estime la significancia estadistica del valor observado en el caso de la red real
'''

#se generan las 16 redes que surjen de colocar 'm' o 'f' a los 'na' de la red original

naNum = np.sum(dolphinsGender[:,1]=='NA')
naGenderPos = list(itertools.product(['m','f'], repeat=naNum))

edgesHeteroRatioNa = []
edgesHeteroRatioRndNa = []
pValNa = []
QNa =[]

for iterNa in range(0,len(naGenderPos)):
    dolphinsNa = dolphins.copy()
    na0=-1
    for n in dolphinsNa.nodes():
        nGender = dolphinsGender[dolphinsGender[:,0]==n,1][0]
        if nGender=='NA':
            na0+=1
            nGender = naGenderPos[iterNa][na0]
        dolphinsNa.nodes[n]["gender"] = nGender
   
    edgesHetero=0
    
    for (e0,e1) in dolphinsNa.edges():
        if dolphinsNa.node[e0]['gender'] != dolphinsNa.node[e1]['gender']:
            edgesHetero+=1
     
    #aca calculo la fraccion de enlaces entre generos diferentes para la red de estudio
    edgesHeteroRatio=edgesHetero/edgesTot
    
    #guardo para cada una de las redes la fraccion de enlaces entre generos distintos
    edgesHeteroRatioNa.append(edgesHeteroRatio)
    
    #estrategia 1: permutacion 

    #luego hago 10000 asignaciones al azar del género a la red dolphinsrnd
    dolphinsRnd=dolphinsNa.copy()
    dolphinsGenderRnd = np.array(dolphinsRnd.nodes("gender"))
    edgesHeteroRatioRnd=[]
    numiter=10000
   
    for i in range(numiter):
        dolphinsGenderRnd[:,1] = random.sample(list(dolphinsGenderRnd[:,1]), nodesTot)
        for n in dolphinsRnd.nodes():
            dolphinsRnd.nodes[n]["gender"] = dolphinsGenderRnd[dolphinsGenderRnd[:,0]==n,1][0]
        edgesHeteroRnd=0
        for (e0,e1) in dolphinsRnd.edges():
            if dolphinsRnd.node[e0]['gender'] != dolphinsRnd.node[e1]['gender']:
                edgesHeteroRnd+=1
        edgesHeteroRatioRnd.append(edgesHeteroRnd/edgesTot)

    # Calculamos la probabilidad de que la proporcion de enlaces entre generos 
    # distintos sea tan o mas extrema que la observada para la red de estudio 
    # considerando que HO es verdadera (considerando como distribucion nula a la 
    # generada por permutar los sexos de los delfines, dejando inalterada la 
    # topologia de la red)

    pVal=1-sum(i >= edgesHeteroRatio for i in edgesHeteroRatioRnd)/numiter
    pValNa.append(pVal)
    edgesHeteroRatioRndNa.append(edgesHeteroRatioRnd)
    
    # Estos p-valores estarian indicando que la red tiene una proporcion de enlaces 
    # entre generos distintos mucho menor que lo esperado por azar por ende, es 
    # homofilica. Esto se mantiene asi para todas las posibles asignaciones de 
    # generos de los delfines cuyo genero no fue definido.
    
    #para cada una de las 16 redes se calcula su modularidad
    
    edgesMod = 0
    for (nodoi, nodoj) in list(dolphins.edges()): 
        # Se recorren todos los enlaces de la red.
        ki=dolphins.degree(nodoi)
        kj=dolphins.degree(nodoj)
        edgesMod += (ki*kj)/(2*edgesTot)
    Q = (edgesHetero - edgesMod)/edgesTot
    
    QNa.append(Q)
    #da un valor positivo, es decir que hay mas enlaces entre nodos del mismo tipo que lo
    #que se esperaria por azar 

#convierto listas resultantes en numpy arrays
    
edgesHeteroRatioNa = np.array(edgesHeteroRatioNa)
edgesHeteroRatioRndNa = np.array(edgesHeteroRatioRndNa)
pValNa = np.array(pValNa)
QNa = np.array(QNa)

#print(dolphins.edges("Ripplefluke"))

#%% 
#plot de las figuras RndNa
#graficamos la distribucion nula (random) para la fraccion de enlaces entre 
#generos diferentes

fig, axs = plt.subplots(nrows=4,ncols=4)
plt.tight_layout()
axs = axs.ravel()
for iterNa in range(0,len(naGenderPos)):
    plt.tight_layout()
    plt.sca(axs[iterNa])
    plt.hist(edgesHeteroRatioRndNa[iterNa,:], bins=40)
    plt.axvline(edgesHeteroRatioNa[iterNa], color='k', linestyle='dashed', linewidth=1)
    plt.title('p-value: ' + str(pValNa[iterNa]) + '[' + str(naGenderPos[iterNa]) + ']',fontsize=10)
    plt.xlim(0,1)
    if iterNa==0:
        plt.axvline(color='k', linestyle='dashed', label='Fraccion original')

fig.text(0.5, 0.025, 'Fraccion de enlaces heterofilicos', ha='center', va='center', fontsize=25)
fig.text(0.20, 0.5, 'Frecuencia', ha='center', va='center', rotation='vertical', fontsize=25)
fig.suptitle('Distribuciones nulas de redes con asignacion de generos al azar', ha='center', va='center', fontsize=25)
fig.legend(bbox_to_anchor=(0.85, 0.5), shadow=True, fontsize='xx-large')

plt.show()

#%%

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


proportions_ztest(edgesHetero, edgesTot, true_mu, 'two-sided', true_var)
   
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


#%% Ej. 2(c)
'''
Identifique alguna metodología basada en observables topológicos para eliminar
nodos secuencialmente de la red de manera de dividirla en dos componentes de tamaños
comparables en el menor número de pasos. Explique y muestre los resultados obtenidos.
Intente cuantificar su estrategia comparándola con lo que se obtendría al eliminar nodos
de manera aleatoria.
'''

# La estrategia sera ir eliminando en cada paso aquel nodo que resulte ser el 
# que constituya parte de la mayor cantidad de "caminos mas cortos posibles" 
# entre dos nodos cualesquiera. Dicha iteracion concluira una vez que la red se 
# halle disconexa. 
# Es esperable que dicha estrategia corte rapidamente la red, aunque debe 
# notarse que esto no garantiza generalmente que las componentes disconexas 
# resulten ser dos de tamaños comparables, por ejemplo no ocurriria asi si el 
# grafo fuera una estrella.

nx.transitivity(dolphins)
dolphinsCut = dolphins.copy()

node2extract=[]
while nx.is_connected(dolphinsCut):
    shortPaths = nx.shortest_path(dolphinsCut)
    numShortPaths = dict((n,0) for n in dolphinsCut.nodes())
    for n1 in dolphinsCut.nodes():
        for n2 in dolphinsCut.nodes():
            shortPathList = shortPaths[n1][n2][1:-1]
            for n in dolphinsCut.nodes():
                if n in shortPathList:
                    numShortPaths[n] = numShortPaths[n]+1
    node2extract.append(max(numShortPaths, key=numShortPaths.get))
    dolphinsCut.remove_node(node2extract[-1])

# Ilustramos el resultado con dicho metodo
dolphinNet = [dolphins,dolphinsCut]
fig, axs = plt.subplots(1,2)
axs = axs.ravel()
for sp in range(0,2):
    plt.sca(axs[sp])
    nx.draw(dolphinNet[sp],
            width=1,
            edge_color = 'c',
            node_color=["blue" if g=="m" else ("red" if g=="f" else "yellow") for g in nx.get_node_attributes(dolphinNet[sp], "gender").values()], 
            node_size=80,
            font_size=10,
            with_labels=True,
           )
    male = mpatches.Patch(color='b', label='Male')
    female = mpatches.Patch(color='r', label='Female')
    unknown = mpatches.Patch(color='y', label='Unknown')
    plt.suptitle('Subtracted nodes: ' ': ' + str(node2extract) + ' (' + str(len(node2extract)) + ' steps)')
    if sp==0:
        plt.title('Original Network')
    elif sp==1:
        plt.title('Cut Network')
        plt.legend(handles=[male,female,unknown])
plt.show()

# Ahora, compararemos con el numero de pasos requeridos al extraer nodos 
# aleatorios. Notese que, para conservar la condicion de separar la red en 
# grafos de tamano comparable, dicha aleatoriedad sera restringida por 
# la condicion de que ninguno de los subgrafos resultantes contenga menos de 10
# nodos. Es decir, se rechazara la eliminacion de nodos perifericos que separen
# a la red en un bloque muy grande y otro muy pequeno.

R = 150
numNodes = []
for r in range(0,R):
    dolphinsCut = dolphins.copy()
    node2extract2=[]
    thresh=0
    try:
        while nx.is_connected(dolphinsCut):
            candidate = random.choice(list(dolphinsCut.nodes()))
            dolphinsCutCand = dolphinsCut.copy()
            dolphinsCutCand.remove_node(candidate)
            subGr = list(nx.connected_component_subgraphs(dolphinsCutCand))
            if len(subGr)>1 and len(min(subGr,key=len))<5 and thresh<10:
                thresh+=1
                continue
            else:
                thresh=0
                node2extract2.append(candidate)
                dolphinsCut = dolphinsCutCand.copy()
        numNodes.append(len(node2extract2))
    except:
        continue
# Graficamos el histograma resultante y comparamos con el metodo anterior
plt.figure()
plt.hist(numNodes, bins=20, color='g', edgecolor='k', label='Random node removal')
plt.axvline(4, color='k', linestyle='dashed', linewidth=1, label='Bridge node removal')
plt.axvline(len(dolphins.nodes()), color='c', linestyle='dashed', linewidth=1, label='Total number of nodes')
plt.legend()
plt.xlabel('Number of nodes')
plt.ylabel('Occurence')
plt.show()

# Se observa que el numero de nodos a extraer es siempre superior en el caso
# aleatorio.
#%% Ej. 3
'''
3) Considere la red as-22july06.gml creada por Mark Newman que contiene la
 estructura de los sistemas autónomos de internet relevada a mediados de 2006.
'''

# Cargar datos

Newman = nx.read_gml(path + 'as-22july06.gml')

#%% Ej. 3(a)
'''
a) Encuentre gráficamente la distribución de grado P_k como función de k explorando
diferentes alternativas: un bineado lineal o logarítmico, utilizando escalas logarítmicas o
lineales en uno o ambos ejes. Discuta que alternativa permite apreciar mejor el carácter
libre de escala de dicha distribución.
'''

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

#%%
'''
b. Utilizando funcionalidad de la librería igraph, estime el exponente de dicha distribución.
'''

from scipy.optimize import curve_fit

x=deg
y=cnt

def modelo(x, m, b):
    return b*(1/x)**m

parametros_iniciales = [1,9500]

popt, pcov = curve_fit(modelo, x,y, p0 = parametros_iniciales)

# Resultado del ajuste:
m=popt[0]
b=popt[1]

print('\nResultado del ajuste (y = b*(1/x)^m):\nm = %2.4f\tb = %2.4f'%(m,b))

x_modelo = np.linspace(1,2500, 10000)

plt.figure()
plt.tight_layout()
plt.plot(x, y,'o',label='Datos')
plt.xscale('log')
plt.yscale('log')
plt.title(r'Ajuste $y=b \left( \dfrac{1}{x^{m}} \right)$')
plt.plot(x_modelo, modelo(x_modelo, *popt), 'r-', label=r'Ajuste')
plt.legend(loc='best',fontsize=LegendSize)
plt.ylabel(r'Cantidad de nodos',fontsize=AxisLabelSize)
plt.xlabel(r'Grado',fontsize=AxisLabelSize)
plt.grid()
plt.tight_layout()
plt.show()


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
'''

redesStr = ['netscience','as-22july06']
redes = {}
avnd = {}
degree0 = {}
nodes = {}
degreeAvnd = {}
ds = {}
d = {}
indx = []

for s in range(len(redesStr)): 
    redes[s] = nx.read_gml(path + redesStr[s] + '.gml')
    nodes[s] = redes[s].nodes()
#nombres de los nodos
    avnd[s] = nx.average_neighbor_degree(redes[s]) 
#diccionario con el grado promedio de los vecinos de cada nodo.
    degree0[s] = redes[s].degree() 
#grado de cada nodo.

    ds[s] = [dict(degree0[s]), avnd[s]]
    
    for k in avnd[s].keys():
       d[s][k] = list(d[s][k] for d[s] in ds[s])

    degreeAvnd[s] = d[s].values()
    degreeAvnd[s] = np.array(list(degreeAvnd[s]))
    
#primero voy a eliminar los 0 de degreeAvnd para poder realizar el ajuste lineal de logx-logy en el punto ii.

    for k in range(len(degreeAvnd[s][:,0])):
        if degreeAvnd[s][k,0] == 0:
            indx.append(k)
            
#unicamente la primer red tiene nodos desconectados, borro esos registros de degreeAvnd

degreeAvnd[0] = np.delete(degreeAvnd[0], indx, axis=0)

#%%
'''
ii. Analizar la tendencia observada en un gráfico que consigne dicho valor k nn (k)
como función del grado.
iii.Asumiendo que k_nn (k) = ak^μ , estime el exponente de correlación a partir de
realizar una regresión de log k_nn ~ log k. Asegurese de graficar el fiteo en el
grafico anterior. 
'''

model = {}
redesTitle = ['colaboraciones cientificas', 'internet']
mu = []
b = []

for s in range(len(degreeAvnd)):

    # x from 0 to 30
    
    x = np.transpose(degreeAvnd[s][:,0])
    y = np.transpose(degreeAvnd[s][:,1])
    
    x = x.reshape((-1,1))  # conversion between (N,) dimension arrays and (N,1) 
    y = y.reshape((-1,1))
    
    logx = np.log10(x)
    logy = np.log10(y)
    
    logx = logx.reshape((-1,1))
    logy = logy.reshape((-1,1))
    
    model[s] = LinearRegression()
    model[s].fit(logx, logy)
    
    logx_new = np.linspace(min(logx), max(logx), 100)
    logy_new = model[s].predict(logx_new[:, np.newaxis])
    
    mu.append(round(float(model[s].coef_),4))
    b.append(round(float(model[s].intercept_),4))
    
    plt.figure()
    plt.plot(x, y,'.k')
    plt.plot(10**logx_new, 10**logy_new,'r', label = r'$\mu$'+ ' = ' + str(mu[s]))
    plt.text(0.5, 0.5, str(mu[s]))
    plt.xlabel(r'Degree, ($k$)', fontsize=20)
    plt.ylabel(r'Average degree of neighbors, ($k_{nn}$)', fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Red de ' + redesTitle[s] ,fontsize=30)
    plt.legend(fontsize=20)
    plt.show()

#%%
'''
iv. Considere la red de colaboraciones y la de internet nuevamente Encuentre
cuantitativamente la asortatividad de la red utilizando ahora el estimador
propuesto por Newman (ver PDF).
Para ello tenga encuenta lo desarrollado en las eqs [8.26 – 8.29] del libro de
Newman.Como se corresponde este coeficiente con el estimado en el punto
anterior? A qué se debe?
'''

rNewman_redes = []
rNX_redes = []
redesTitle = ['colaboraciones cientificas', 'internet']

for s in range(len(redes)): 
    # Se arma un vector 'netDeg4' con los grados de cada nodo:
    netDeg4 = np.array(redes[s].degree())
    # Lista con los grados de cada nodo
    netDeg4Grados = netDeg4[:,1].astype(int)   
    
    # Usando el método del libro (Newman, ecuaciones 8.26-8.29):
    
    S1=sum(netDeg4Grados)
    S2=sum(netDeg4Grados**2)
    S3=sum(netDeg4Grados**3)
    
    suma=0
    for (nodoi, nodoj) in list(redes[s].edges()): 
    # Se recorren todos los enlaces de la red.
        ki=redes[s].degree(nodoi)
        kj=redes[s].degree(nodoj)
        suma += ki*kj
    
    Se=2*suma
    
    rNewman=(S1*Se-S2**2)/(S1*S3-S2**2)
    rNewman_redes.append(rNewman)
    
    
    rNX=nx.degree_assortativity_coefficient(redes[s])
    rNX_redes.append(rNX)
    
    print('Red de ' + redesTitle[s])
    print('Coeficiente de correlación (Newman):\nr= ',rNewman)
    print('Coeficiente de correlación (función de Networkx):\nr= ',rNX)

# Falta comparar con el punto anterior.

#%%
'''
b) Corra el script de cálculo (puntos i-iii) para las redes Y2H y AP-MS. Puede explicar lo que observa en cuanto a la asortatividad reportada?
'''
'''
i.Determine, para nodos de grado k, cuánto vale en media el grado de sus vecinos. 
'''

redesStr = ['Y2H','AP-MS']
redes = {}
avnd = {}
degree0 = {}
nodes = {}
degreeAvnd = {}
ds = {}
d = {}
indx = []

for s in range(len(redesStr)):
    redes[s] = nx.Graph(ldata(path + 'yeast_' + redesStr[s] + '.txt'))
    nodes[s] = redes[s].nodes()
    avnd[s] = nx.average_neighbor_degree(redes[s]) 
    degree0[s] = redes[s].degree() 
    
    ds[s] = [dict(degree0[s]), avnd[s]]
    
    for k in avnd[s].keys():
       d[s][k] = list(d[s][k] for d[s] in ds[s])

    degreeAvnd[s] = d[s].values()
    degreeAvnd[s] = np.array(list(degreeAvnd[s]))

#%%
'''
ii. Analizar la tendencia observada en un gráfico que consigne dicho valor k nn (k)
como función del grado.
iii.Asumiendo que k_nn (k) = ak^μ , estime el exponente de correlación a partir de
realizar una regresión de log k_nn ~ log k. Asegurese de graficar el fiteo en el
grafico anterior. 
'''

model = {}
redesTitle = ['Yeast Two-Hybrid', 'Affinity Purification - Mass Spectrometry']
mu = []
b = []

for s in range(len(degreeAvnd)):

    # x from 0 to 30
    
    x = np.transpose(degreeAvnd[s][:,0])
    y = np.transpose(degreeAvnd[s][:,1])
    
    x = x.reshape((-1,1))  # conversion between (N,) dimension arrays and (N,1) 
    y = y.reshape((-1,1))
    
    logx = np.log10(x)
    logy = np.log10(y)
    
    logx = logx.reshape((-1,1))
    logy = logy.reshape((-1,1))
    
    model[s] = LinearRegression()
    model[s].fit(logx, logy)
    
    logx_new = np.linspace(min(logx), max(logx), 100)
    logy_new = model[s].predict(logx_new[:, np.newaxis])
    
    mu.append(round(float(model[s].coef_),4))
    b.append(round(float(model[s].intercept_),4))
    
    plt.figure()
    plt.plot(x, y,'.k')
    plt.plot(10**logx_new, 10**logy_new,'r', label = r'$\mu$'+ ' = ' + str(mu[s]))
    plt.text(0.5, 0.5, str(mu[s]))
    plt.xlabel(r'Degree, ($k$)', fontsize=20)
    plt.ylabel(r'Average degree of neighbors, ($k_{nn}$)', fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Red de ' + redesTitle[s] ,fontsize=30)
    plt.legend(fontsize=20)
    plt.show()
    
#%%
