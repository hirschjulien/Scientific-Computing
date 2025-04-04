import numpy as np
import matplotlib.pyplot as plt


"""
Ns: number of realisations
T: Time steps
N: number of intervals
"""

def WienerProcess(T, N, Ns, seed):
    np.random.seed(seed)
    dt = T / N
    dW = np.sqrt(dt) * np.random.randn(Ns, N)
    W = np.hstack([np.zeros((Ns, 1)), np.cumsum(dW, axis=1)])
    Tw = np.linspace(0, T, N+1)

    return W, Tw, dW

Ns = 10
W, Tw, dW = WienerProcess(10, 1000, Ns, 100)
# plt.plot(Tw, W[1])
print(Tw)
print(np.shape(W))

def euler_expexp(ffun, gfun, T, x0, W, args):
    N = len(T) - 1
    nx = len(x0)
    X = np.zeros((nx, N+1))

    X[:,0] = x0

    for i in range(1,N):
        dt = T[i+1] - T[i]
        dW = W[:,i+1] - W[:,i]
        f = ffun(T[i], X[:,i], args[0])
        g = gfun(T[i], X[:,i], args[1])
        X[:, i+1] = X[:,i] + f * dt + g * dW

    return X

def vdP(t, x, mu):
    xdot = np.array([x[1], mu * (1 - x[0]**2) * x[1] - x[0]])
    return xdot

def vdP_g(t, x, sigma):
    return np.array([0.0, sigma * (1.0 + x[0] * x[0])])

x0 = np.array([0.5, 0.5])
mu = 3
sigma = 1

X = np.zeros(np.shape(W))
for i in range(Ns):
    X[i,:] = euler_expexp(vdP, vdP_g, Tw, x0, W, [mu, sigma])