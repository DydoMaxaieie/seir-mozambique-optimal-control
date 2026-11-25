"""Capítulo 4: agrupamento de Louvain da matriz de fluxos entre províncias."""
import numpy as np, networkx as nx
from collections import Counter
T=np.array([
[0,986,119,74,24,13,17,16,12,46,23],
[396,0,1547,509,172,68,145,85,62,267,229],
[111,1513,0,399,66,32,77,45,39,139,113],
[196,214,1476,0,150,127,401,56,58,235,180],
[61,73,89,310,0,626,434,30,21,81,43],
[2,0,1,45,196,0,624,67,15,62,21],
[16,1,54,515,66,646,0,97,42,111,55],
[0,3,1,3,2,0,2,0,868,710,115],
[0,2,3,1,0,2,1,897,0,779,129],
[76,184,172,223,117,82,122,128,148,0,12447],
[838,2275,1797,2471,1307,985,1400,1753,2006,21166,0]],float)
nm=['Niassa','Cabo Delgado','Nampula','Zambézia','Tete','Manica','Sofala','Inhambane','Gaza','Maputo','Maputo Cidade']
W=T+T.T
G=nx.Graph()
for i in range(11):
    for j in range(i+1,11): G.add_edge(i,j,weight=W[i,j])
def lab(c): return tuple(sorted(tuple(sorted(nm[i] for i in s)) for s in c))
for res in [0.6,0.8,1.0,1.2,1.5,2.0]:
    cnt=Counter(); Q=[]
    for seed in range(200):
        c=nx.community.louvain_communities(G,weight='weight',resolution=res,seed=seed)
        cnt[lab(c)]+=1; Q.append(nx.community.modularity(G,c,weight='weight',resolution=res))
    print('res',res,'n particoes',len(cnt),'Qmedio',round(np.mean(Q),3))
    for p,n in cnt.most_common(3): print('  ',n,'/200',[list(g) for g in p])
print('--- variantes a resolucao 1')
def run(Wm,directed=False,tag=''):
    G=nx.DiGraph() if directed else nx.Graph()
    for i in range(11):
        for j in range(11):
            if i!=j and (directed or i<j): G.add_edge(i,j,weight=Wm[i,j])
    cnt=Counter()
    for seed in range(200):
        c=nx.community.louvain_communities(G,weight='weight',seed=seed); cnt[lab(c)]+=1
    print(tag,len(cnt))
    for p,n in cnt.most_common(2): print('  ',n,[list(g) for g in p])
P=T/T.sum(1,keepdims=True)
run(P+P.T,False,'prob. linha simetrizada')
run(np.maximum(P,P.T),False,'max(P,P^T)')
run(T,True,'dirigido fluxos')
run(P,True,'dirigido prob.')
run(np.sqrt(np.outer(T.sum(1),T.sum(0)))*0+T/np.sqrt(np.outer(T.sum(1)+T.sum(0),T.sum(1)+T.sum(0))),False,'normalizado por grau')
