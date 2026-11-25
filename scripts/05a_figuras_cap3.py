import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from seirmoz import verif_escalar as lq, verif_seir as seir
from seirmoz.paths import RESULTS
col=['#1B365D','#4E7D5B','#B07A1F','#8C2F39','#6B5B95']
# ---- LQ ----
t,x,q,u,it=lq.fbs('hermite',80)
tt=np.linspace(0,1,400)
fig,ax=plt.subplots(1,3,figsize=(12,3.6))
ax[0].plot(tt,lq.xe(tt),color=col[0],lw=1.8,label='Solução exacta $x(t)$')
ax[0].plot(t[::8],x[::8],'o',color=col[3],ms=4,label='FBS com RK4')
ax[0].plot(tt,0.5*np.exp(tt),color='#9AA5B1',ls='--',lw=1.2,label='Sem controlo')
ax[0].set_xlabel('$t$',fontsize=9);ax[0].set_title('Estado',fontsize=10);ax[0].legend(fontsize=8,frameon=False)
ax[1].plot(tt,lq.qe(tt),color=col[0],lw=1.8,label='$q(t)$ exacto')
ax[1].plot(t[::8],q[::8],'o',color=col[3],ms=4,label='$q(t)$ FBS')
ax[1].plot(tt,-lq.qe(tt),color=col[1],lw=1.4,ls='-.',label='$u^*(t)=-q(t)$')
ax[1].axhline(0,color='#9AA5B1',lw=0.8);ax[1].set_xlabel('$t$',fontsize=9);ax[1].set_title('Adjunta e controlo',fontsize=10);ax[1].legend(fontsize=8,frameon=False)
hs=1/np.array([10,20,40,80,160]);lab={'euler':'Euler','rk4avg':'RK4, ponto médio por média','hermite':'RK4, ponto médio por Hermite'}
E={}
for m,c in zip(('euler','rk4avg','hermite'),col):
    e=[np.abs(lq.fbs(m,n)[1]-lq.xe(np.linspace(0,1,n+1))).max() for n in (10,20,40,80,160)];E[m]=e
    ax[2].loglog(hs,e,'o-',color=c,lw=1.4,ms=4,label=lab[m])
for p,c in ((1,col[0]),(2,col[1]),(4,col[2])): pass
ax[2].set_xlabel('Passo $h$',fontsize=9);ax[2].set_ylabel('Erro máximo em $x$',fontsize=9);ax[2].set_title('Ordem de convergência',fontsize=10);ax[2].legend(fontsize=7.5,frameon=False)
for a in ax: a.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig(str(RESULTS/'verificacao'/'exemplo_escalar.png'),dpi=170);plt.close()
for m in E: print(m,['%.2e'%v for v in E[m]])
# ---- SEIR ----
T=seir.T;tt=np.arange(seir.n+1)*seir.h
Bs=[0.25,0.5,1,2];res={}
X0,J0=seir.fwd(np.zeros(seir.n+1),1.0,1.0)
fig,ax=plt.subplots(2,2,figsize=(10.5,6.4))
ax[0,0].plot(tt,X0[:,2],color='#8C2F39',lw=1.8,label='Sem controlo')
for B,c in zip(Bs,[col[0],col[1],col[2],col[4]]):
    u,X,Q,J,it,hist=seir.solve(1.0,B,w=0.3,tol=1e-7,it=800);res[B]=(u,X,Q,J,it)
    ax[0,0].plot(tt,X[:,2],color=c,lw=1.4,label=f'$B={B}$')
    ax[0,1].plot(tt,u,color=c,lw=1.4,label=f'$B={B}$')
ax[0,0].set_title('Infecciosos $I(t)$ (fracção da população)',fontsize=10);ax[0,1].set_title('Controlo óptimo $u^*(t)$',fontsize=10)
ax[0,0].legend(fontsize=8,frameon=False);ax[0,1].legend(fontsize=8,frameon=False)
u,X,Q,J,it=res[1]
for j,(nm,c) in enumerate(zip(['$q_S$','$q_E$','$q_I$'],col)): ax[1,0].plot(tt,Q[:,j],color=c,lw=1.5,label=nm)
ax[1,0].set_title('Variáveis adjuntas ($B=1$)',fontsize=10);ax[1,0].legend(fontsize=8,frameon=False)
for w,c in ((0.5,col[3]),(0.3,col[0])):
    uu,XX,QQ,JJ,ii,hist=seir.solve(1.0,0.25,w=w,tol=1e-7,it=120)
    ax[1,1].semilogy(np.arange(1,len(hist)+1),hist,color=c,lw=1.5,label=f'$\\varpi={w}$')
ax[1,1].set_title('Convergência do FBS ($B=0.25$)',fontsize=10);ax[1,1].set_xlabel('Iteração',fontsize=9);ax[1,1].set_ylabel('Variação máxima do controlo',fontsize=9);ax[1,1].legend(fontsize=8,frameon=False)
for a in ax[:,:].ravel(): a.spines[['top','right']].set_visible(False)
for a in ax[0]: a.set_xlabel('Dias',fontsize=9)
ax[1,0].set_xlabel('Dias',fontsize=9)
fig.tight_layout();fig.savefig(str(RESULTS/'verificacao'/'exemplo_seir.png'),dpi=170);plt.close()
print('J0',J0,'pico0',X0[:,2].max(),'t',tt[X0[:,2].argmax()],'ataque0',1-X0[-1,0])
for B in Bs:
    u,X,Q,J,it=res[B];print(B,'pico',X[:,2].max().round(4),'t',tt[X[:,2].argmax()],'ataque',(1-X[-1,0]).round(3),'J/J0',(J/J0).round(3),'umed',u.mean().round(3),'umax',u.max().round(3),'it',it)
# saturacao com caixa [-1,1] no LQ
t,x,q,u,it=lq.fbs('hermite',80,lim=1.0); print('LQ caixa [-1,1]: u=-1 ate t =',t[np.where(np.abs(u[::2])>=0.9999)[0].max()] if (np.abs(u[::2])>=0.9999).any() else None, 'J',lq.J_of(x,u[::2],t))
t,x,q,u,it=lq.fbs('hermite',80); print('J sem saturacao',lq.J_of(x,u[::2],t),'q(0)',q[0],'u(0)',u[0],'iters',it)
