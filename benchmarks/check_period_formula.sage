import sys; sys.path.insert(0, os.getcwd())
from schottky_theta import *
K=ComplexField(100); i=K(0,1)
D={1:Disc(-2-2*i,K(1)),-1:Disc(-2+2*i,K(1)),2:Disc(2-2*i,K(1)),-2:Disc(2+2*i,K(1))}
S=SchottkyTheta(D,40); F=S.solve(); P=S.period_matrix(F)
print(P.apply_map(lambda x:x.n(digits=12)))
def cr(a,b,c,d): return (a-c)*(b-d)/((a-d)*(b-c))
def mu(n):
    A,B=S.fixed[n]; a,b,c,d=S.sigma[n].list(); return (a*d-b*c)/(c*B+d)**2
for L in [5,7,9,10]:
    B=Matrix(K,2,2)
    for m in (1,2):
        for n in (1,2):
            tot=K(0)
            def rec(g,wd):
                global tot
                if wd and abs(wd[0])!=m and abs(wd[-1])!=n:
                    Am,Bm=S.fixed[m]; An,Bn=S.fixed[n]
                    tot+=cr(Bm,Am,mobius(g,Bn),mobius(g,An)).log()
                if len(wd)==L: return
                for j in S.I:
                    if wd and wd[-1]==-j: continue
                    rec(g*S.sigma[j],wd+[j])
            rec(S.sigma[1]**0,[])
            B[m-1,n-1]=tot+(mu(n).log() if m==n else cr(S.fixed[m][1],S.fixed[m][0],S.fixed[n][1],S.fixed[n][0]).log())
    print(L,(B-P).apply_map(abs).apply_map(lambda x:x.n(digits=3)))
