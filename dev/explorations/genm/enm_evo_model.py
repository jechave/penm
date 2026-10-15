import numpy as np
from scipy.optimize import minimize
X=np.load('ca.npy'); N=len(X)
D=np.linalg.norm(X[:,None]-X[None,:],axis=2)
I,J=np.where(np.triu(D<14.0,1)); M=len(I); cidx=np.arange(M); d0=D[I,J].copy()
K0,RC,W=0.25,10.5,1.0                     # k fitted to the ddG scale
def kof(l): return 0.5*K0*(1.0-np.tanh((l-RC)/W))
SIG=0.87; A_STATES=20                     # sigma fitted to structural divergence
rs=np.random.default_rng(2024)
U0=rs.normal(0,SIG,(A_STATES,M)); U1=rs.normal(0,SIG,(A_STATES,M))
s0=np.zeros(N,int)
def lengths(seq): return d0+(U0[seq[I],cidx]-U0[s0[I],cidx])+(U1[seq[J],cidx]-U1[s0[J],cidx])
def V(x,L,k):
    x=x.reshape(N,3); dd=np.linalg.norm(x[J]-x[I],axis=1); return 0.5*np.sum(k*(dd-L)**2)
def G(x,L,k):
    x=x.reshape(N,3); r=x[J]-x[I]; dd=np.linalg.norm(r,axis=1); e=r/dd[:,None]
    c=(k*(dd-L))[:,None]*e; g=np.zeros((N,3)); np.add.at(g,I,-c); np.add.at(g,J,c); return g.ravel()
def mini(x0,L,k):
    r=minimize(V,x0,args=(L,k),jac=G,method='L-BFGS-B',options={'maxiter':50000,'ftol':1e-18,'gtol':1e-12})
    return r.x,r.fun
def hess(x,L,k):
    x=x.reshape(N,3); r=x[J]-x[I]; dd=np.linalg.norm(r,axis=1); e=r/dd[:,None]; g=(dd-L)/dd
    P=e[:,:,None]*e[:,None,:]; B=-k[:,None,None]*(P+g[:,None,None]*(np.eye(3)-P))
    H=np.zeros((N,3,N,3))
    np.add.at(H,(I,slice(None),J,slice(None)),B); np.add.at(H,(J,slice(None),I,slice(None)),B)
    np.add.at(H,(I,slice(None),I,slice(None)),-B); np.add.at(H,(J,slice(None),J,slice(None)),-B)
    return H.reshape(3*N,3*N)
def kabsch(Y,Z):
    Y=Y.reshape(N,3); Z=Z.reshape(N,3); cy=Y.mean(0); cz=Z.mean(0)
    U,S,Vt=np.linalg.svd((Y-cy).T@(Z-cz)); ds=np.sign(np.linalg.det(Vt.T@U.T))
    return (((Vt.T@np.diag([1,1,ds])@U.T)@(Y-cy).T).T+cz).flatten()
