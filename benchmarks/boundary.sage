# Sup-norm error of Omega_1 on the boundary of the fundamental domain (64 points on each of the 4 circles)
# and at the interior point 0.3+0.1i, for the genus 2 datum of Section 7 (discs of radius 1 at +-2+-2i).
import sys, time, os
sys.path.insert(0, os.getcwd())
from schottky_theta import *
def config(K):
    i = K(0,1); r = K(1)
    return {1:Disc(-2-2*i,r),-1:Disc(-2+2*i,r),2:Disc(2-2*i,r),-2:Disc(2+2*i,r)}
def points(K):
    pts = []
    for D in config(K).values():
        for k in range(64):
            th = 2*K.pi()*k/64
            pts.append(D.center + D.radius*K(th.cos(), th.sin()))
    return pts
Kr = ComplexField(500)
Sref = SchottkyTheta(config(Kr), 150); Fr = Sref.solve()
zin = Kr(3)/10 + Kr(0,1)/10
ref_b = [Sref.Omega(1, Fr.column(0), z) for z in points(Kr)]
ref_i = Sref.Omega(1, Fr.column(0), zin)
print('reference done', flush=True)
for p in [10, 20, 40, 60, 80]:
    K = ComplexField(int(0.4*p*3.33) + 60)
    S = SchottkyTheta(config(K), p); F = S.solve()
    eb = max(abs(Kr(S.Omega(1, F.column(0), z)) - r) for z, r in zip(points(K), ref_b))
    ei = abs(Kr(S.Omega(1, F.column(0), K(3)/10 + K(0,1)/10)) - ref_i)
    print('p=%d  sup error on the boundary circles %.1e   error at 0.3+0.1i %.1e' % (p, float(eb), float(ei)), flush=True)
