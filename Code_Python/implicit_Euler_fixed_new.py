import numpy as np 
import matplotlib.pyplot as plt
#from euler_explicit import *

# xinit is the first guess, xk is the step before that
def Newton(fun, jac, tk, xk, dt, xinit, tol, maxit, *args):
    k = 0
    t = tk + dt
    x = xinit
    f = fun(t, x, *args)
    J = jac(t, x, *args)
    R = x - f * dt - xk
    I = np.eye(len(xk))
    while(k < maxit and np.linalg.norm(R) > tol):
        k = k + 1
        dRdx = I - J * dt
        dx = np.linalg.solve(dRdx, R)
        x = x - dx
        f = fun(t, x, *args)
        J = jac(t, x, *args)
        R = x - dt * f - xk
    
    return x


def euler_im(fun, jac, t0, T, N, x0, *args):
    dt = (T - t0) / N
    xn = len(x0)
    y = np.zeros((xn,N+1))
    t = np.zeros(N+1)

    tol = 1.0e-8
    maxit = 100

    t[0] = t0
    y[:,0] = x0

    for i in range(N):
        f = fun(t[i], y[:,i], *args)
        t[i+1] = t[i] + dt
        yinit = y[:,i] + f * dt #explicit Euler as intial guess
        y[:, i+1] = Newton(fun, jac, t[i], y[:,i], dt, yinit, tol, maxit, *args)
    
    return t, y

def vdP(t, x, mu):
    xdot = np.zeros(2)
    xdot[0] = x[1]
    xdot[1] = mu * (1 - x[0]*x[0]) * x[1] - x[0]
    return xdot

def jacvdP(t, x, mu):
    jac = np.zeros((2,2))
    jac[0,1] = 1.0
    jac[1,0] = -2* mu * x[0] * x[1] - 1.0
    jac[1,1] = mu * (1 - x[0] * x[0])
    return jac

x0 = [2.0, 0.0]
mu = 1
N = 2000
t0 = 0
T = 50

#sol_ex = euler_ex(vdP, t0, T, N, x0, mu)
sol_im = euler_im(vdP, jacvdP, t0, T, N, x0, mu)
#plt.plot(sol_ex[1][0], sol_ex[1][1])
#plt.show()
plt.plot(sol_im[1][0], sol_im[1][1])
#plt.title(f'Implicit, N={N}')