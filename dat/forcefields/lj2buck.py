from math import exp


def lj2buck(epsilon, sigma, delta=12.):
    a = 6*epsilon*exp(delta)/(delta-6)
    rho = sigma*pow(2, 1/6.)/delta
    c = epsilon*delta*pow(sigma*pow(2, 1/6.), 6.)/(delta-6)
    return a, rho, c


def add_buckingham(f):
    for pt in f.particle_types:
        a, rho, c = lj2buck(pt.epsilon, pt.sigma)
        pt.a = a
        pt.rho = rho
        pt.c = c

    return f
