def load_ref(Kr):
    while not os.path.exists('benchmarks/server/ref.txt'): time.sleep(20)
    re,im=open('benchmarks/server/ref.txt').read().split()
    R=RealBallField(500)
    return Kr(R(QQ(re)),R(QQ(im)))
