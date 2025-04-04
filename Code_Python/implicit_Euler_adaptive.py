import numpy as np
import matplotlib.pyplot as plt

def vdP(t, x, mu):
    xdot = np.array([x[1], mu * (1 - x[0]**2) * x[1] - x[0]])
    return xdot

def jacvdP(t, x, mu):
    jac = np.array([[0, 1], [-2 * mu * x[0] * x[1] - 1, mu * (1 - x[0]**2)]])
    return jac

def NewtonODE(fun, jac, tk, xk, dt, xinit, maxit, mu): 
    k = 0 #start iterating index, counts iterations
    t = tk + dt #the time one step after the initial time 
    x = xinit
    f = fun(t,x,mu) #evaluate f and Jacobian at initial guess x0
    J = jac(t,x,mu)
    R = x - dt*f - xk #actual value minus predicted value
    I = np.eye(len(xk))
    while k<maxit:
        k = k + 1
        dRdx = I - J*dt
        dx = np.linalg.solve(dRdx, R)
        x = x - dx
        f = fun(t,x,mu) #evaluate f and Jacobian at initial guess x0
        J = jac(t,x,mu)
        R = x - dt*f - xk
    return x

def ImplicitEulerAdaptiveStepSize(fun, jac, ta, tb, dt, x0, mu):
    eps = 0.8
    abstol = 1e-3
    reltol = 1e-3
    maxit = 50

    T = [ta]
    X = [x0]
    t = ta

    while t < tb:
        if t + dt > tb:
            dt = tb - t

        xk = X[-1]  # current x is the last element in X (x_k = X[-1])
        f = fun(t, xk, mu)
        J = jac(t, xk, mu)
        xinit = xk + f * dt  

        # Full step
        X1 = NewtonODE(fun, jac, t, xk, dt, xinit, maxit, mu)

        # Two half steps
        dt_half = dt / 2  
        xinit_half = xk + dt_half * f
        X_half = NewtonODE(fun, jac, t, xk, dt_half, xinit_half, maxit, mu)
        X2 = NewtonODE(fun, jac, t + dt_half, X_half, dt_half, xinit_half, maxit, mu)

        e = X2 - X1
        r = np.max(np.abs(e) / np.maximum(abstol, np.abs(X2) * reltol))

        if r <= 1.0:
            X.append(X2)
            t += dt
            T.append(t)

        dt = max(1e-6, min(2.0 * dt, np.sqrt(eps / r) * dt))  

    return np.array(T), np.array(X).T

sol = ImplicitEulerAdaptiveStepSize(vdP, jacvdP, 0, 50, 1, np.array([2, 0]), 10)
plt.plot(sol[0], sol[1][0])
plt.plot(sol[0], sol[1][1])
plt.show()




# import numpy as np
# import matplotlib.pyplot as plt

# def vdP(t, x, mu):
#     xdot = np.zeros(2)
#     xdot[0] = x[1]
#     xdot[1] = mu * (1 - x[0]*x[0]) * x[1] - x[0]
#     return xdot

# def jacvdP(t, x, mu):
#     jac = np.zeros((2,2))
#     jac[0,1] = 1.0
#     jac[1,0] = -2* mu * x[0] * x[1] - 1.0
#     jac[1,1] = mu * (1 - x[0] * x[0])
#     return jac

# def NewtonODE(fun, jac, tk, xk, dt, xinit, maxit, mu): 
#     k = 0 #start iterating index, counts iterations
#     t = tk + dt #the time one step after the initial time 
#     x = xinit
#     f = fun(t,x,mu) #evaluate f and Jacobian at initial guess x0
#     J = jac(t,x,mu)
#     R = x - dt*f - xk #actual value minus predicted value
#     I = np.eye(len(xk))
#     while k<maxit:
#         k = k + 1
#         dRdx = I - J*dt
#         dx = np.linalg.solve(dRdx, R)
#         x = x - dx
#         f = fun(t,x,mu) #evaluate f and Jacobian at initial guess x0
#         J = jac(t,x,mu)
#         R = x - dt*f - xk
#     return x

# def ImplicitEulerAdaptiveStepSize(fun, jac, ta, tb, dt, x0, mu):
#     eps = 0.8
#     #dt = (tb - ta) / N

#     ###### no inital N, start with inital step size #####
#     #nx = len(x0)
#     #X = np.zeros((nx, N+1))
#     #T = np.zeros(N+1)

#     abstol = 1.0e-2
#     reltol = 1.0e-2
#     maxit = 50

#     # Eulers Implicit Method
#     #T[0] = ta
#     #X = x0
#     #t = ta

#     T = np.array([ta])
#     X = np.array([x0])
#     t = ta

#     while t < tb:
#         if t + dt > tb:
#             dt = tb - t

#         # xk = X[-1]
#         f = fun(t, X[-1], mu)
#         J = jac(t, X[-1], mu)
#         xinit = X[-1] + f * dt

#         # Full step
#         X1 = NewtonODE(fun, jac, t, X[-1], dt, xinit, maxit, mu)

#         # Two half steps
#         hm = dt / 2
#         dt_half = t + hm
#         xinit_half = X[-1] + hm * f
#         X_half = NewtonODE(fun, jac, t, X[-1], dt_half, xinit_half, maxit, mu)
#         X2 = NewtonODE(fun, jac, t + dt_half, X_half, dt_half, xinit_half, maxit, mu)

#         # Error estimation
#         e = X2 - X1
#         r = np.max(np.abs(e) / np.maximum(abstol, np.abs(X2) * reltol))
        
#         if r <= 1.0:
#             X = np.vstack([X, X2])
#             #T[t+1] = T[t] + dt
#             t += dt
#             T = np.append(T, t)
#         else:
#             dt = max(1e-6, min(2.0 * dt, np.sqrt(eps / r) * dt))

#     return T, X



# sol = ImplicitEulerAdaptiveStepSize(vdP, jacvdP, 0, 1, 5, np.array([1,0]), 10)
# print(sol[1])
# plt.plot(sol[1][0], sol[1][1])




