import autograd.numpy as np
from autograd import grad
from math import exp

def f(x1, x2):
    y = x1**2*np.exp(x2) + np.cos(x2)**2
    return y
fGrad = grad(f,(0,1))

print(fGrad(3.,np.pi))