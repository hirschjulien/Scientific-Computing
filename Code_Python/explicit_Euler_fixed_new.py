import numpy as np 
import matplotlib.pyplot as plt

def euler_ex(fun, t0, T, N, x0, *args):
    dt = (T - t0) / N
    xn = len(x0)
    y = np.zeros((xn,N+1))
    t = np.zeros(N+1)

    t[0] = t0
    y[:,0] = x0

    for i in range(N):
        f = fun(t[i], y[:,i], *args)
        t[i+1] = t[i] + dt
        y[:,i+1] = y[:,i] + f * dt
    
    return t, y


# van der Pol
def vdP(t, x, mu):
    xdot = np.zeros(2)
    xdot[0] = x[1]
    xdot[1] = mu * (1 - x[0]*x[0]) * x[1] - x[0]
    return xdot

x0 = [2.0, 0.0]
mu = 3
N = 5000
t0 = 0
T = 10

sol = euler_ex(vdP, 0, 10, N, x0, mu)
plt.plot(sol[0], sol[1][1])
plt.title(f'Explicit, N={N}')