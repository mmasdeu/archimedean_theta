import sys,time,os
sys.path.insert(0, os.getcwd())   # run from the repository root
from schottky_theta import *
def config(K,r=1):
    i=K(0,1); r=K(r)
    return {1:Disc(-2-2*i,r),-1:Disc(-2+2*i,r),2:Disc(2-2*i,r),-2:Disc(2+2*i,r)}
