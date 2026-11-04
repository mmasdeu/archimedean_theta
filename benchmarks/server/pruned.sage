# usage: sage pruned.sage d exponent   (mu = 10^-exponent)
load('benchmarks/server/common.sage')
d=int(sys.argv[-2]); e=int(sys.argv[-1])
Kr=ComplexBallField(500)
load('benchmarks/server/refload.sage')
K=ComplexField(int(d*3.33)+40); z=K(3)/10+K(0,1)/10
S=SchottkyTheta(config(K),3)
t0=time.time(); a,c=adaptive_Omega(S,1,z,K(10)^(-e)); tt=time.time()-t0
print('pruned d=%d mu=1e-%d terms=%d time=%.1f (summation done)'%(d,e,c,tt),flush=True)
ref=load_ref(Kr)
print('pruned d=%d mu=1e-%d terms=%d time=%.1f err=%.1e'%(d,e,c,tt,float(abs(Kr(a)-ref))),flush=True)
