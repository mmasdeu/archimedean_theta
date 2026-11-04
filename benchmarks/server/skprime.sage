import sys,time,os
sys.path.insert(0, os.getcwd())   # run from the repository root
from theta_arch import *
CC=ComplexField(400); i=CC.gen()
balls={1:Ball(CC(1)/2,CC(1)/10),2:Ball(CC.gen()/2,CC(1)/10)}
z=(-1-i)/2; w=CC(1)
for prec in [10,20,40,60]:
    t0=time.time(); X=SKprime(balls,z,w,prec); print('SKprime p=%d time=%.2f X=%s'%(prec,time.time()-t0,X.n(digits=30)),flush=True)
