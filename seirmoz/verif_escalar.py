"""Problema linear-quadrático escalar com solução exacta (Secção 3.3.1)."""
import numpy as np
s3=np.sqrt(3)
c1=0.5*(1+s3)*np.exp(-s3)/((1+s3)*np.exp(-s3)-(1-s3)*np.exp(s3)); c2=0.5-c1
xe=lambda t:c1*np.exp(s3*t)+c2*np.exp(-s3*t)
qe=lambda t:(1-s3)*c1*np.exp(s3*t)+(1+s3)*c2*np.exp(-s3*t)
f=lambda x,u:x+u
qd=lambda x,q:-(2*x+q)
def J_of(x,u,t): return np.trapezoid(x**2+0.5*u**2,t)
def fbs(method,n,tol=1e-12,it=800,w=0.5,lim=2.0):
    h=1.0/n; t=np.linspace(0,1,n+1)
    if method=='hermite':
        th=np.linspace(0,1,2*n+1); u=np.zeros(2*n+1)
    else: u=np.zeros(n+1)
    for k in range(it):
        x=np.empty(n+1);x[0]=.5
        if method=='euler':
            for i in range(n): x[i+1]=x[i]+h*f(x[i],u[i])
        elif method=='rk4avg':
            for i in range(n):
                um=(u[i]+u[i+1])/2; xi=x[i]
                k1=f(xi,u[i]);k2=f(xi+h/2*k1,um);k3=f(xi+h/2*k2,um);k4=f(xi+h*k3,u[i+1]); x[i+1]=xi+h/6*(k1+2*k2+2*k3+k4)
        else:
            for i in range(n):
                ua,um,ub=u[2*i],u[2*i+1],u[2*i+2]; xi=x[i]
                k1=f(xi,ua);k2=f(xi+h/2*k1,um);k3=f(xi+h/2*k2,um);k4=f(xi+h*k3,ub); x[i+1]=xi+h/6*(k1+2*k2+2*k3+k4)
        q=np.empty(n+1);q[n]=0
        if method=='euler':
            for i in range(n,0,-1): q[i-1]=q[i]-h*qd(x[i],q[i])
            un=np.clip(-q,-lim,lim)
        elif method=='rk4avg':
            for i in range(n,0,-1):
                xm=(x[i]+x[i-1])/2;um=(u[i]+u[i-1])/2;qi=q[i]
                k1=qd(x[i],qi);k2=qd(xm,qi-h/2*k1);k3=qd(xm,qi-h/2*k2);k4=qd(x[i-1],qi-h*k3); q[i-1]=qi-h/6*(k1+2*k2+2*k3+k4)
            un=np.clip(-q,-lim,lim)
        else:
            fx=np.array([f(x[i],u[2*i]) for i in range(n+1)]); xm=(x[:-1]+x[1:])/2+h/8*(fx[:-1]-fx[1:])
            for i in range(n,0,-1):
                # adjunta com u e x conhecidos: q' depende so de x,q
                qi=q[i]
                k1=qd(x[i],qi);k2=qd(xm[i-1],qi-h/2*k1);k3=qd(xm[i-1],qi-h/2*k2);k4=qd(x[i-1],qi-h*k3); q[i-1]=qi-h/6*(k1+2*k2+2*k3+k4)
            gq=np.array([qd(x[i],q[i]) for i in range(n+1)]); qm=(q[:-1]+q[1:])/2+h/8*(gq[:-1]-gq[1:])
            un=np.empty(2*n+1); un[0::2]=np.clip(-q,-lim,lim); un[1::2]=np.clip(-qm,-lim,lim)
        d=np.abs(un-u).max(); u=w*un+(1-w)*u
        if d<tol: break
    return t,x,q,u,k+1
if __name__=='__main__':
    for m in ('euler','rk4avg','hermite'):
        errs=[]
        for n in (10,20,40,80,160):
            t,x,q,u,it=fbs(m,n); errs.append(np.abs(x-xe(t)).max())
        ord_=[np.log2(errs[i]/errs[i+1]) for i in range(len(errs)-1)]
        print(m,['%.2e'%e for e in errs],np.round(ord_,2),'iters',it)
    t,x,q,u,it=fbs('hermite',80); print('q(0)',q[0],qe(0),'u max',np.abs(u).max())
