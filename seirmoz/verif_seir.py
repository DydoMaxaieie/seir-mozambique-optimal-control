"""Distanciamento óptimo no SEIR homogéneo (Secção 3.3.2)."""
import numpy as np
sg,gm,beta=1/5.2,1/10,0.2; T=180.0; umax=0.6; h=0.25; n=int(T/h)
x0=np.array([1-1e-3,5e-4,5e-4,0.0])
def f(x,u):
    S,E,I,R=x; lam=beta*(1-u)*I
    return np.array([-lam*S,lam*S-sg*E,sg*E-gm*I,gm*I])
def adj(q,x,u,A):
    S,E,I,R=x; lam=beta*(1-u)*I; qS,qE,qI,qR=q
    return np.array([lam*(qS-qE), sg*(qE-qI), -A-beta*(1-u)*S*(qE-qS)+gm*(qI-qR), 0.0])
def fwd(u,A,B):
    X=np.empty((n+1,4));X[0]=x0
    for i in range(n):
        x=X[i];ua,ub=u[i],u[i+1];um=(ua+ub)/2
        k1=f(x,ua);k2=f(x+h/2*k1,um);k3=f(x+h/2*k2,um);k4=f(x+h*k3,ub);X[i+1]=x+h/6*(k1+2*k2+2*k3+k4)
    J=np.trapezoid(A*X[:,2]+0.5*B*u**2,dx=h);return X,J
def bwd(X,u,A):
    Q=np.zeros((n+1,4))
    for i in range(n,0,-1):
        q=Q[i];xa,xb=X[i-1],X[i];xm=(xa+xb)/2;um=(u[i]+u[i-1])/2
        k1=adj(q,xb,u[i],A);k2=adj(q-h/2*k1,xm,um,A);k3=adj(q-h/2*k2,xm,um,A);k4=adj(q-h*k3,xa,u[i-1],A);Q[i-1]=q-h/6*(k1+2*k2+2*k3+k4)
    return Q
def solve(A,B,w=0.5,tol=1e-8,it=500,u0=None):
    u=np.zeros(n+1) if u0 is None else u0.copy();hist=[]
    for k in range(it):
        X,J=fwd(u,A,B);Q=bwd(X,u,A)
        S,E,I,R=X.T;un=np.clip((Q[:,1]-Q[:,0])*beta*S*I/B,0,umax)
        d=np.abs(un-u).max();hist.append(d);u=w*un+(1-w)*u
        if d<tol: break
    X,J=fwd(u,A,B);Q=bwd(X,u,A);return u,X,Q,J,k+1,hist
def main():

    t=np.arange(n+1)*h
    X0,J0=fwd(np.zeros(n+1),1.0,1.0)
    print('sem controlo: pico I',X0[:,2].max().round(4),'dia',t[X0[:,2].argmax()],'final',(1-X0[-1,0]).round(3),'J',round(J0,3))
    for B in (0.1,0.5,2,8):
        u,X,Q,J,it,hist=solve(1.0,B)
        print(B,'pico',X[:,2].max().round(4),'ataque',(1-X[-1,0]).round(3),'J',round(J,3),'J/J0',round(J/J0,3),'umedio',u.mean().round(3),'umax',u.max().round(3),'it',it)
    # convergencia com relaxacao
    for w in (1.0,0.7,0.5,0.3):
        u,X,Q,J,it,hist=solve(1.0,8.0,w=w,it=300); print('w',w,'it',it,'d final',hist[-1])
    # gradiente: adjunta vs dif. finitas
    u=np.full(n+1,0.2);X,J=fwd(u,1.0,8.0);Q=bwd(X,u,1.0);S,E,I,R=X.T
    gr=h*(8.0*u-(Q[:,1]-Q[:,0])*beta*S*I)
    for i in (60,200,400):
        e=1e-6;up=u.copy();up[i]+=e;um=u.copy();um[i]-=e
        fd=(fwd(up,1.0,8.0)[1]-fwd(um,1.0,8.0)[1])/(2*e);print(i,'FD',fd,'adj',gr[i],'rel',abs(fd-gr[i])/abs(fd))


if __name__=='__main__':
    main()
