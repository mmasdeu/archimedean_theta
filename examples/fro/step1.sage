load('examples/fro/fro.sage')
K = ComplexField(200)
eta, theta = K.pi()/12, K.pi()/4
for g in [2,3]:
    discs, gens = fro_datum(eta, theta, g, K)
    check_orientation(discs, gens)
