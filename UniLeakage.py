#! /usr/bin/env python
# Arunas L Radzvilavicius
# ucbprad >at< ucl.ac.uk
# GNU GPL v3

from numpy import *
#from matplotlib.pyplot import *
#from matplotlib.mlab import griddata
import sys
from scipy.special import binom
def Run(LKx,lkx,mux,xix):
    IM=-1; t=xix; M=50; mu=mux;
    LK = LKx;
    lk = lkx;
    m = arange(M+1)[:,None]; w = 1.0 - (1.0*m/M)**t
    U = array( [ [ binom(M-j,i-j)*(mu)**(i-j)*(1-mu)**(M-i) for j in range(M+1)] for i in range(M+1) ] )
    S = array( [ [ 1.0*binom(2*j,i)*binom(2*M-2*j,M-i)/binom(2*M,M) for j in range(M+1)] for i in range(M+1) ] )
    F1 = array( [ [ 1.0*binom(2*j,i)*binom(4*M-2*j,2*M-i)/binom(4*M,2*M) for j in range(2*M+1)] for i in range(2*M+1) ] )
    F2 = array( [ [ 1.0*binom(j,i)*binom(2*M-j,M-i)/binom(2*M,M) for j in range(2*M+1)] for i in range(2*M+1) ] )
    FF = dot(F2,F1)
    m1 = int(LK*M)
    #Q1 = array( [ [ 1.0*binom(j,i)*binom(M-j,m1-i)/binom(M,m1) for j in range(M+1)] for i in range(m1+1) ] )
    Q1 = array( [ [ binom(m1,i)*(1.0*j/M)**(i)*(1-1.0*j/M)**(m1-i) for j in range(M+1)] for i in range(m1+1) ] )	
    m0 = 2*M-int(LK*M)
    Q0 = array( [ [ binom(m0,i)*(1.0*j/M)**(i)*(1-1.0*j/M)**(m0-i) for j in range(M+1)] for i in range(m0+1) ] )
    #Q0=array( [ [ 1.0*binom(j,i-j)*binom(M-j,m0-M-i+j)/binom(M,m0-M) if i>=j else 0 for j in range(M+1)] for i in range(m0+1) ] )	
    m1 = int(lk*M)
    #q1=array( [ [ 1.0*binom(j,i)*binom(M-j,m1-i)/binom(M,m1) for j in range(M+1)] for i in range(m1+1) ] )
    q1 = array( [ [ binom(m1,i)*(1.0*j/M)**(i)*(1-1.0*j/M)**(m1-i) for j in range(M+1)] for i in range(m1+1) ] ) 
    m0 = 2*M-int(lk*M)
    q0 = array( [ [ binom(m0,i)*(1.0*j/M)**(i)*(1-1.0*j/M)**(m0-i) for j in range(M+1)] for i in range(m0+1) ] )
    #q0=array( [ [ 1.0*binom(j,i-j)*binom(M-j,m0-M-i+j)/binom(M,m0-M) if i>=j else 0 for j in range(M+1)] for i in range(m0+1) ] ) 
    def Omega(s):
        t = dot( diagflat(w), s ) 
        b = dot( w.T, sum(s, axis=1) )[0]
        return 1.0*t/b
    s0 = zeros((M+1,2)) # MatingType1
    s1 = zeros((M+1,2)) # MatingType2
    s0[M/2,0] = 1
    s1[M/2,0] = 1
    p01 = [];
    maxgen = int(1e8)
    for i in range(maxgen):
        if i==200:
            f = 0.01
            ds = f*s0[:,0]
            s0[:,0] -= ds 
            s0[:,1] += ds
        p01.append( sum(s0[:,1]) )
	s0 = dot(U,s0); s1 = dot(U,s1);
        s0 = Omega(s0); s1 = Omega(s1);
        g0 = 1.0*s0; g1 = 1.0*s1;
        z00 = convolve( dot(Q0,g0[:,0]),dot(Q1,g1[:,0]) )[:,None]
        z01 = convolve( dot(Q0,g0[:,0]),dot(Q1,g1[:,1]) )[:,None]
        z10 = convolve( dot(q0,g0[:,1]),dot(q1,g1[:,0]) )[:,None]
        z11 = convolve( dot(q0,g0[:,1]),dot(q1,g1[:,1]) )[:,None]
        s0[:,0] = dot(FF,0.5*z00+0.5*z01)[0:M+1].flatten()
        s0[:,1] = (0.5*dot(FF,z11)[0:M+1]+0.5*dot(FF,z10)[0:M+1]).flatten()
        s1[:,0] = (0.5*dot(FF,z00)[0:M+1]+0.5*dot(FF,z10)[0:M+1]).flatten()
        s1[:,1] = (0.5*dot(FF,z11)[0:M+1]+0.5*dot(FF,z01)[0:M+1]).flatten()
        s0 *= 2; s1 *= 2;
	if i>1500 and p01[-1]>0.99999: break
	if i>1500 and p01[-1]<0.00001: break	
        if i>1500 and abs(p01[-10]-p01[-1])<1e-8: break
        
    return p01[-1]

LK = float(sys.argv[1])
lk = float(sys.argv[2])
mu = float(sys.argv[3])
xi = float(sys.argv[4])
if LK>1 or lk>1: sys.exit(0)
print LK, lk, Run(LK,lk,mu,xi)
