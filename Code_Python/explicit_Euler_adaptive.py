import numpy as np 
import matplotlib.pyplot as plt

def ExplicitEulerAdaptiveStepSize(fun, ta, tb, dt, x0, *args):
    eps = 0.8
    X = [x0]
    T = [ta]
    t = ta

    abstol = 1e-5
    reltol = 1e-5

    naccept = 0
    nreject = 0

    while t < tb:
        if t + dt > tb:
            dt = tb - t

        xk = X[-1]
        f = fun(t, xk, *args)

        X1 = xk + dt * f

        hm = dt / 2
        tm = t + hm
        xm = xk + hm * f
        fm = fun(tm, xm, *args)
        X2 = xm + hm * fm

        e = X2 - X1
        r = np.max(np.abs(e) / np.maximum(abstol, np.abs(X2) * reltol))

        if r <= 1.0:
            X.append(X2)
            t += dt
            T.append(t)
            naccept += 1

        dt = max(1e-6, min(2.0 * dt, np.sqrt(eps / r) * dt))
        nreject += 1   

    return np.array(T), np.array(X).T, naccept, nreject


# van der Pol
def vdP(t, x, mu):
    xdot = np.zeros(2)
    xdot[0] = x[1]
    xdot[1] = mu * (1 - x[0]*x[0]) * x[1] - x[0]
    return xdot

x0 = [2.0, 0.0]
mu = 10
ta = 0
tb = 50
dt = 5

sol = ExplicitEulerAdaptiveStepSize(vdP, ta, tb, dt, x0, mu)
plt.plot(sol[0], sol[1][0])
plt.plot(sol[0], sol[1][1])
plt.show

print(f"naccept: {sol[2]} \nnreject: {sol[3]}")