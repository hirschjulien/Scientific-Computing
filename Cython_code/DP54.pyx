import cython
import numpy as np
cimport numpy as cnp
from libc.stdlib cimport malloc, free
from libc.math cimport fabs, pow, fmax, fmin, exp
from libc.stdio cimport printf
import time 
from libc.string cimport memcpy

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void fill_array(double* X, double[:] x0, int n_cells) noexcept nogil:
    cdef int i, j
    for i in range(n_cells):
        for j in range(3):
            X[i*3 + j] = x0[j]

ctypedef (double, double, double) (*function_type)(double[3], double[:]) noexcept nogil
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (double, double, double) functionvan(double[3] x, double[:] argv) noexcept nogil:
    cdef double[3] f
    f[0] = x[1]
    f[1] = argv[0] * (1 - x[0] * x[0]) * x[1] - x[0]
    f[2] = 0.0
    return f[0], f[1], f[2]
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (double, double, double) functionprey(double[3] x, double[:] argv) noexcept nogil:
    cdef double[3] f
    f[0] = argv[0] * (1 - x[1]) * x[0]
    f[1] = - argv[1] * (1 - x[0]) * x[1]
    f[2] = 0.0
    return f[0], f[1], f[2]
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (double, double, double) functionCSTR3(double[3] x, double[:] argv) noexcept nogil:
    cdef double[3] f
    cdef double k = argv[0] * exp(- argv[1] / x[2])    # k0 * exp(-Ea_R / T)
    cdef double r = k * x[0] * x[1]
    cdef double beta = - argv[2] / (argv[3] * argv[4])
    cdef double FdV = (argv[5] / argv[6]) 
    f[0] = (argv[7] - x[0]) * FdV - r       # (CA_in - CA) * (F / V) - r  
    f[1] = (argv[8] - x[1]) * FdV - 2 * r   # (CB_in - CB) * (F / V) - 2.0 * r 
    f[2] = (argv[9] - x[2]) * FdV + r * beta  #(T_in - T) * (F / V) + beta * r
    return f[0], f[1], f[2]
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (double, double, double) functionCSTR1(double[3] x, double[:] argv) noexcept nogil:
    cdef double[3] f
    cdef double k = argv[0] * exp(- argv[1] / x[2])    # k0 * exp(-Ea_R / T)
    cdef double r = k * argv[7] * argv[8]
    cdef double beta = - argv[2] / (argv[3] * argv[4])
    cdef double FdV = (argv[5] / argv[6]) 
    f[0] = 0
    f[1] = 0
    f[2] = (argv[9] - x[2]) * FdV + r * beta  #(T_in - T) * (F / V) + beta * r
    return f[0], f[1], f[2]
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (double, double, double) functionstability(double[3] y, double[:] argv) noexcept nogil:
    cdef double[3] f
    # y[0]Re，y[1]Im
    # argv[0]λRe，argv[1]λIm    
    # (a+bi)(c+di) = (ac-bd) + (ad+bc)i
    f[0] = argv[0] * y[0] - argv[1] * y[1]  
    f[1] = argv[0] * y[1] + argv[1] * y[0]  
    f[2] = 0
    return f[0], f[1], f[2]


ctypedef (double, double, double) (*function_type2)(double[3], double[3], double[3], double[:]) noexcept nogil  
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (double, double, double) functionPFR3(double[3] x, double[3] x_l, double[3] x_r, double[:] argv) noexcept nogil:
    cdef double[3] f
    cdef double k0 = argv[0]      
    cdef double Ea_R = argv[1]    
    cdef double deltaHr = argv[2] 
    cdef double rho = argv[3]     
    cdef double cP = argv[4]      
    cdef double v = argv[5]       
    cdef double D_A = argv[6]     
    cdef double D_B = argv[7]     
    cdef double D_T = argv[8]     

    cdef double dz   = argv[10]
    cdef double k = k0 * exp(-Ea_R / x[2])
    cdef double r = k * x[0] * x[1]
    cdef double beta = -deltaHr / (rho * cP)
    cdef double dz2 = dz * dz

    cdef double dCA_dz = (x[0] - x_l[0]) / dz
    cdef double dCB_dz = (x[1] - x_l[1]) / dz
    cdef double dT_dz  = (x[2] - x_l[2]) / dz

    cdef double d2CA_dz2 = (x_l[0] - 2.0 * x[0] + x_r[0]) / (dz2)
    cdef double d2CB_dz2 = (x_l[1] - 2.0 * x[1] + x_r[1]) / (dz2)
    cdef double d2T_dz2  = (x_l[2]  - 2.0 * x[2] + x_r[2] ) / (dz2)

    f[0] = - v * dCA_dz + D_A * d2CA_dz2 - r
    f[1] = - v * dCB_dz + D_B * d2CB_dz2 - 2.0 * r
    f[2]  = - v * dT_dz  + D_T * d2T_dz2  + beta * r

    return f[0], f[1], f[2]
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline (double, double, double) functionPFR1(double[3] x, double[3] x_l, double[3] x_r, double[:] argv) noexcept nogil:
    cdef double[3] f
    cdef double k0 = argv[0]      
    cdef double Ea_R = argv[1]    
    cdef double deltaHr = argv[2] 
    cdef double rho = argv[3]     
    cdef double cP = argv[4]      
    cdef double v = argv[5]       
    cdef double D_A = argv[6]     
    cdef double D_B = argv[7]     
    cdef double D_T = argv[8]     

    cdef double dz   = argv[10]
    cdef double k = k0 * exp(-Ea_R / x[2])
    cdef double r = k * argv[11] * argv[12]
    cdef double beta = -deltaHr / (rho * cP)
    cdef double dz2 = dz * dz

    cdef double dT_dz  = (x[2] - x_l[2]) / dz

    cdef double d2T_dz2  = (x_l[2]  - 2.0 * x[2] + x_r[2] ) / (dz2)

    f[0] = 0
    f[1] = 0
    f[2]  = - v * dT_dz  + D_T * d2T_dz2  + beta * r

    return f[0], f[1], f[2]

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void DP54_step_PFR(double* X, double[:] argv, double h, double* x_new, double* x_hat, double[21] a, double[21] k, double[7] b, double[7] b_hat, function_type2 func) noexcept nogil:
    cdef double[3] temp
    cdef int N  = <int>argv[9]
    cdef int n = 0
    cdef double[3] x, x_r
    cdef double[3] x_l
    x_l[0] = argv[11]
    x_l[1] = argv[12]
    x_l[2] = argv[13]
    for n in range(N):
        if n == N-1:
            x_r[0] = argv[14]   # CA_out
            x_r[1] = argv[15]   # CB_out
            x_r[2] = argv[16]   # T_out
        else:
            x_r[0] = X[(n+1) * 3]
            x_r[1] = X[(n+1) * 3 + 1]
            x_r[2] = X[(n+1) * 3 + 2]
        x[0] = X[n * 3]
        x[1] = X[n * 3 + 1]
        x[2] = X[n * 3 + 2]
        k[0], k[1], k[2] = func(x, x_l, x_r, argv)
        temp[0] = x[0] + h * a[0] * k[0]
        temp[1] = x[1] + h * a[0] * k[1]
        temp[2] = x[2] + h * a[0] * k[2]

        k[3], k[4], k[5] = func(temp, x_l, x_r, argv)
        temp[0] = x[0] + h * a[1] * k[0] + h * a[2] * k[3]
        temp[1] = x[1] + h * a[1] * k[1] + h * a[2] * k[4]
        temp[2] = x[2] + h * a[1] * k[2] + h * a[2] * k[5]

        k[6], k[7], k[8] = func(temp, x_l, x_r, argv)
        temp[0] = x[0] + h * a[3] * k[0] + h * a[4] * k[3] + h * a[5] * k[6]
        temp[1] = x[1] + h * a[3] * k[1] + h * a[4] * k[4] + h * a[5] * k[7]
        temp[2] = x[2] + h * a[3] * k[2] + h * a[4] * k[5] + h * a[5] * k[8]

        k[9], k[10], k[11] = func(temp, x_l, x_r, argv)
        temp[0] = x[0] + h * a[6] * k[0] + h * a[7] * k[3] + h * a[8] * k[6] + h * a[9] * k[9]
        temp[1] = x[1] + h * a[6] * k[1] + h * a[7] * k[4] + h * a[8] * k[7] + h * a[9] * k[10]
        temp[2] = x[2] + h * a[6] * k[2] + h * a[7] * k[5] + h * a[8] * k[8] + h * a[9] * k[11]

        k[12], k[13], k[14] = func(temp, x_l, x_r, argv)
        temp[0] = x[0] + h * a[10] * k[0] + h * a[11] * k[3] + h * a[12] * k[6] + h * a[13] * k[9] + h * a[14] * k[12]  
        temp[1] = x[1] + h * a[10] * k[1] + h * a[11] * k[4] + h * a[12] * k[7] + h * a[13] * k[10] + h * a[14] * k[13]  
        temp[2] = x[2] + h * a[10] * k[2] + h * a[11] * k[5] + h * a[12] * k[8] + h * a[13] * k[11] + h * a[14] * k[14] 

        k[15], k[16], k[17] = func(temp, x_l, x_r, argv)
        temp[0] = x[0] + h * a[15] * k[0] + h * a[16] * k[3] + h * a[17] * k[6] + h * a[18] * k[9] + h * a[19] * k[12] + h * a[20] * k[15]  
        temp[1] = x[1] + h * a[15] * k[1] + h * a[16] * k[4] + h * a[17] * k[7] + h * a[18] * k[10] + h * a[19] * k[13] + h * a[20] * k[16]  
        temp[2] = x[2] + h * a[15] * k[2] + h * a[16] * k[5] + h * a[17] * k[8] + h * a[18] * k[11] + h * a[19] * k[14] + h * a[20] * k[17]

        k[18], k[19], k[20] = func(temp, x_l, x_r, argv)

        x_new[n * 3 ] = x[0] + h * (b[0] * k[0] + b[1] * k[3] + b[2] * k[6] + b[3] * k[9] + b[4] * k[12] + b[5] * k[15] + b[6] * k[18])
        x_new[n * 3 + 1] = x[1] + h * (b[0] * k[1] + b[1] * k[4] + b[2] * k[7] + b[3] * k[10] + b[4] * k[13] + b[5] * k[16] + b[6] * k[19])
        x_new[n * 3 + 2] = x[2] + h * (b[0] * k[2] + b[1] * k[5] + b[2] * k[8] + b[3] * k[11] + b[4] * k[14] + b[5] * k[17] + b[6] * k[20])

        x_hat[n * 3 ] = x[0] + h * (b_hat[0] * k[0] + b_hat[1] * k[3] + b_hat[2] * k[6] + b_hat[3] * k[9] + b_hat[4] * k[12] + b_hat[5] * k[15] + b_hat[6] * k[18])
        x_hat[n * 3 + 1] = x[1] + h * (b_hat[0] * k[1] + b_hat[1] * k[4] + b_hat[2] * k[7] + b_hat[3] * k[10] + b_hat[4] * k[13] + b_hat[5] * k[16] + b_hat[6] * k[19])
        x_hat[n * 3 + 2] = x[2] + h * (b_hat[0] * k[2] + b_hat[1] * k[5] + b_hat[2] * k[8] + b_hat[3] * k[11] + b_hat[4] * k[14] + b_hat[5] * k[17] + b_hat[6] * k[20])
        x_l[0] = X[3 * n    ]      
        x_l[1] = X[3 * n + 1]
        x_l[2] = X[3 * n + 2]


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline void DP54_step(double[3] x, double[:] argv, double h, double[3] x_new, double[3] x_hat, double[21] a, double[21] k, double[7] b, double[7] b_hat, function_type func) noexcept nogil:
    cdef double[3] temp
    k[0], k[1], k[2] = func(x, argv)
    temp[0] = x[0] + h * a[0] * k[0]
    temp[1] = x[1] + h * a[0] * k[1]
    temp[2] = x[2] + h * a[0] * k[2]

    k[3], k[4], k[5] = func(temp, argv)
    temp[0] = x[0] + h * a[1] * k[0] + h * a[2] * k[3]
    temp[1] = x[1] + h * a[1] * k[1] + h * a[2] * k[4]
    temp[2] = x[2] + h * a[1] * k[2] + h * a[2] * k[5]

    k[6], k[7], k[8] = func(temp, argv)
    temp[0] = x[0] + h * a[3] * k[0] + h * a[4] * k[3] + h * a[5] * k[6]
    temp[1] = x[1] + h * a[3] * k[1] + h * a[4] * k[4] + h * a[5] * k[7]
    temp[2] = x[2] + h * a[3] * k[2] + h * a[4] * k[5] + h * a[5] * k[8]

    k[9], k[10], k[11] = func(temp, argv)
    temp[0] = x[0] + h * a[6] * k[0] + h * a[7] * k[3] + h * a[8] * k[6] + h * a[9] * k[9]
    temp[1] = x[1] + h * a[6] * k[1] + h * a[7] * k[4] + h * a[8] * k[7] + h * a[9] * k[10]
    temp[2] = x[2] + h * a[6] * k[2] + h * a[7] * k[5] + h * a[8] * k[8] + h * a[9] * k[11]

    k[12], k[13], k[14] = func(temp, argv)
    temp[0] = x[0] + h * a[10] * k[0] + h * a[11] * k[3] + h * a[12] * k[6] + h * a[13] * k[9] + h * a[14] * k[12]  
    temp[1] = x[1] + h * a[10] * k[1] + h * a[11] * k[4] + h * a[12] * k[7] + h * a[13] * k[10] + h * a[14] * k[13]  
    temp[2] = x[2] + h * a[10] * k[2] + h * a[11] * k[5] + h * a[12] * k[8] + h * a[13] * k[11] + h * a[14] * k[14] 

    k[15], k[16], k[17] = func(temp, argv)
    temp[0] = x[0] + h * a[15] * k[0] + h * a[16] * k[3] + h * a[17] * k[6] + h * a[18] * k[9] + h * a[19] * k[12] + h * a[20] * k[15]  
    temp[1] = x[1] + h * a[15] * k[1] + h * a[16] * k[4] + h * a[17] * k[7] + h * a[18] * k[10] + h * a[19] * k[13] + h * a[20] * k[16]  
    temp[2] = x[2] + h * a[15] * k[2] + h * a[16] * k[5] + h * a[17] * k[8] + h * a[18] * k[11] + h * a[19] * k[14] + h * a[20] * k[17]

    k[18], k[19], k[20] = func(temp, argv)

    x_new[0] = x[0] + h * (b[0] * k[0] + b[1] * k[3] + b[2] * k[6] + b[3] * k[9] + b[4] * k[12] + b[5] * k[15] + b[6] * k[18])
    x_new[1] = x[1] + h * (b[0] * k[1] + b[1] * k[4] + b[2] * k[7] + b[3] * k[10] + b[4] * k[13] + b[5] * k[16] + b[6] * k[19])
    x_new[2] = x[2] + h * (b[0] * k[2] + b[1] * k[5] + b[2] * k[8] + b[3] * k[11] + b[4] * k[14] + b[5] * k[17] + b[6] * k[20])

    x_hat[0] = x[0] + h * (b_hat[0] * k[0] + b_hat[1] * k[3] + b_hat[2] * k[6] + b_hat[3] * k[9] + b_hat[4] * k[12] + b_hat[5] * k[15] + b_hat[6] * k[18])
    x_hat[1] = x[1] + h * (b_hat[0] * k[1] + b_hat[1] * k[4] + b_hat[2] * k[7] + b_hat[3] * k[10] + b_hat[4] * k[13] + b_hat[5] * k[16] + b_hat[6] * k[19])
    x_hat[2] = x[2] + h * (b_hat[0] * k[2] + b_hat[1] * k[5] + b_hat[2] * k[8] + b_hat[3] * k[11] + b_hat[4] * k[14] + b_hat[5] * k[17] + b_hat[6] * k[20])

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void compute_PFR(double* t_history, double* X, double* X_history, double[:] argv, double END, double[:] h_size, double reps, double aeps, double* t_adaptive_size, double* t_adaptive, int* numbers_ptr, function_type2 func) noexcept nogil:
    cdef double[7] c = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0]
    cdef double[7] b_hat = [5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, - 92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0]
    cdef double[7] b = [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0]
    cdef double[21] a = [
        # Row 2 (index 0)
        1.0 / 5.0,
        # Row 3 (indices 1, 2)
        3.0 / 40.0, 9.0 / 40.0,
        # Row 4 (indices 3, 4, 5)
        44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0,
        # Row 5 (indices 6, 7, 8, 9)
        19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0,
        # Row 6 (indices 10, 11, 12, 13, 14)
        9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0,
        # Row 7 (indices 15-20)
        35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0
    ]
    cdef int step_counts = 1
    cdef double t_start = 0
    cdef double t_current = t_start
    cdef double h_max = h_size[1] 
    cdef double h_min = h_size[0]
    cdef int i = 0
    cdef double t_end = END
    #cdef double[3] x_new
    #cdef double[3] x_hat
    cdef int total_size = <int>argv[9] * 3
    cdef double *x_new = <double*>malloc(total_size * sizeof(double))
    cdef double *x_hat = <double*>malloc(total_size * sizeof(double))
    cdef double h
    cdef double h_
    cdef double[21] k
    #cdef double[3] err
    #cdef double[3] tol
    cdef double *err = <double*>malloc(total_size * sizeof(double))
    cdef double *tol = <double*>malloc(total_size * sizeof(double))
    cdef double safety = 0.9
    cdef double err_ratio
    cdef int num = 0
    numbers_ptr[0] = 0
    h = h_max 
    while t_current < t_end: 
        if t_current + h > t_end:
            h = t_end - t_current
        DP54_step_PFR(X, argv, h, x_new, x_hat, a, k, b, b_hat, func)
        err[0] = fabs(x_new[0] - x_hat[0])
        tol[0] = reps * fabs(x_new[0]) + aeps
        err_ratio = 0.0
        for num in range(total_size):
            err[num] = fabs(x_new[num] - x_hat[num])
            tol[num] = reps * fabs(x_new[num]) + aeps
            err_ratio = fmax(err[num]/tol[num], err_ratio)
        #err_ratio = fmax(err[0]/tol[0], err[1]/tol[1])
        #err_ratio = fmax(err_ratio, err[2]/tol[2])
        #printf("err_ratio: %f\n", err_ratio)
        if err_ratio <= 1.0:
            t_current = t_current + h
            num = 0
            for num in range(total_size):
                X[num] = x_new[num]
            t_adaptive_size[numbers_ptr[0]] = h
            t_adaptive[numbers_ptr[0]] = t_current
            t_history[numbers_ptr[0]] = t_current
            #x_history_1[numbers_ptr[0]] = x[0]
            #x_history_2[numbers_ptr[0]] = x[1]
            #x_history_3[numbers_ptr[0]] = x[2]
            for num in range(total_size):
                X_history[total_size * numbers_ptr[0] + num] = X[num]
            h *= 2
            #h = h_max
            numbers_ptr[0] += 1
        else:
            h = fmax(safety * h * pow(err_ratio, -0.2), h_min)
            if h <= h_min:
                h = h_min
                t_current += h
                num = 0
                for num in range(total_size):
                    X[num] = x_new[num]
                t_adaptive_size[numbers_ptr[0]] = h
                t_adaptive[numbers_ptr[0]] = t_current
                t_history[numbers_ptr[0]] = t_current
                #x_history_1[numbers_ptr[0]] = x[0]
                #x_history_2[numbers_ptr[0]] = x[1]
                #x_history_3[numbers_ptr[0]] = x[2]
                for num in range(total_size):
                    X_history[total_size * numbers_ptr[0] + num] = X[num]
                numbers_ptr[0] += 1


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void compute(double* t_history, double* x_history_1, double* x_history_2, double* x_history_3, double[3] x, double[:] argv, double END, double[:] h_size, double reps, double aeps, double* t_adaptive_size, double* t_adaptive, int* numbers_ptr, function_type func) noexcept nogil:
    cdef double[7] c = [0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0]
    cdef double[7] b_hat = [5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, - 92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0]
    cdef double[7] b = [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0]
    cdef double[21] a = [
        # Row 2 (index 0)
        1.0 / 5.0,
        # Row 3 (indices 1, 2)
        3.0 / 40.0, 9.0 / 40.0,
        # Row 4 (indices 3, 4, 5)
        44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0,
        # Row 5 (indices 6, 7, 8, 9)
        19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0,
        # Row 6 (indices 10, 11, 12, 13, 14)
        9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0,
        # Row 7 (indices 15-20)
        35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0
    ]
    cdef int step_counts = 1
    cdef double t_start = 0
    cdef double t_current = t_start
    cdef double h_max = h_size[1] 
    cdef double h_min = h_size[0]
    cdef int i = 0
    cdef double t_end = END
    cdef double[3] x_new
    cdef double[3] x_hat
    cdef double h
    cdef double h_
    cdef double[21] k
    cdef double[3] err
    cdef double[3] tol
    cdef double safety = 0.9
    cdef double err_ratio
    numbers_ptr[0] = 0
    h = h_max 
    x[0] = x_history_1[i] 
    x[1] = x_history_2[i]
    x[2] = x_history_3[i]
    while t_current < t_end: 
        if t_current + h > t_end:
            h = t_end - t_current
        DP54_step(x, argv, h, x_new, x_hat, a, k, b, b_hat, func)
        err[0] = fabs(x_new[0] - x_hat[0])
        err[1] = fabs(x_new[1] - x_hat[1])
        err[2] = fabs(x_new[2] - x_hat[2])
        tol[0] = reps * fabs(x_new[0]) + aeps
        tol[1] = reps * fabs(x_new[1]) + aeps
        tol[2] = reps * fabs(x_new[2]) + aeps
        err_ratio = fmax(err[0]/tol[0], err[1]/tol[1])
        err_ratio = fmax(err_ratio, err[2]/tol[2])
        if err_ratio <= 1.0:
            t_current = t_current + h
            x[0] = x_new[0]
            x[1] = x_new[1]
            x[2] = x_new[2]
            t_adaptive_size[numbers_ptr[0]] = h
            t_adaptive[numbers_ptr[0]] = t_current
            t_history[numbers_ptr[0]] = t_current
            x_history_1[numbers_ptr[0]] = x[0]
            x_history_2[numbers_ptr[0]] = x[1]
            x_history_3[numbers_ptr[0]] = x[2]
            h *= 2
            numbers_ptr[0] += 1
        else:
            h = fmax(safety * h * pow(err_ratio, -0.2), h_min)
            if h <= h_min:
                h = h_min
                t_current += h
                x[0] = x_new[0]
                x[1] = x_new[1]
                x[2] = x_new[2]
                t_adaptive_size[numbers_ptr[0]] = h
                t_adaptive[numbers_ptr[0]] = t_current
                t_history[numbers_ptr[0]] = t_current
                x_history_1[numbers_ptr[0]] = x[0]
                x_history_2[numbers_ptr[0]] = x[1]
                x_history_3[numbers_ptr[0]] = x[2]
                numbers_ptr[0] += 1

def DP54_solver(double[:] x0, double[:] argv, double[:] t_range, double[:] h_size, double reps, double aeps, int model, int MAX_STEPS):
    cdef double start = time.perf_counter()
    cdef double[3] x
    cdef int N 
    if model == 5 or model == 6:
        N = <int>argv[9] * 3 
    else:
        N = 1
    cdef int numbers
    cdef double END = t_range[1]
    cdef cnp.ndarray[double, ndim=1] t_history = np.zeros((MAX_STEPS), dtype=np.float64) 
    cdef cnp.ndarray[double, ndim=1] x_history_1 = np.zeros((MAX_STEPS), dtype=np.float64)
    cdef cnp.ndarray[double, ndim=1] x_history_2 = np.zeros((MAX_STEPS), dtype=np.float64)
    cdef cnp.ndarray[double, ndim=1] x_history_3 = np.zeros((MAX_STEPS), dtype=np.float64)
    cdef cnp.ndarray[double, ndim=1] X = np.zeros((N), dtype=np.float64)
    cdef cnp.ndarray[double, ndim=1] X_history = np.zeros((N * MAX_STEPS), dtype=np.float64)
    
    cdef double* t_history_ptr = <double*>t_history.data
    cdef double* x_history_1_ptr = <double*>x_history_1.data
    cdef double* x_history_2_ptr = <double*>x_history_2.data
    cdef double* x_history_3_ptr = <double*>x_history_3.data
    cdef double* X_ptr = <double*>X.data
    cdef double* X_history_ptr = <double*>X_history.data
    cdef cnp.ndarray[double, ndim=1] t_adaptive_size = np.zeros((MAX_STEPS), dtype=np.float64) 
    cdef cnp.ndarray[double, ndim=1] t_adaptive = np.zeros((MAX_STEPS), dtype=np.float64) 
    cdef double* t_adaptive_size_ptr = <double*>t_adaptive_size.data
    cdef double* t_adaptive_ptr = <double*>t_adaptive.data
    
    x_history_1[0] = x0[0]
    x_history_2[0] = x0[1]
    x_history_3[0] = x0[2]
    if model == 1:
        compute(t_history_ptr, x_history_1_ptr, x_history_2_ptr, x_history_3_ptr, x, argv, END, h_size, reps, aeps, t_adaptive_size_ptr, t_adaptive_ptr, &numbers, functionvan)
    elif model == 2:
        compute(t_history_ptr, x_history_1_ptr, x_history_2_ptr, x_history_3_ptr, x, argv, END, h_size, reps, aeps, t_adaptive_size_ptr, t_adaptive_ptr, &numbers, functionprey)
    elif model == 3:
        compute(t_history_ptr, x_history_1_ptr, x_history_2_ptr, x_history_3_ptr, x, argv, END, h_size, reps, aeps, t_adaptive_size_ptr, t_adaptive_ptr, &numbers, functionCSTR3)
    elif model == 4:
        compute(t_history_ptr, x_history_1_ptr, x_history_2_ptr, x_history_3_ptr, x, argv, END, h_size, reps, aeps, t_adaptive_size_ptr, t_adaptive_ptr, &numbers, functionCSTR1)
    elif model == 5:
        fill_array(X_ptr, x0, <int>argv[9])
        compute_PFR(t_history_ptr, X_ptr, X_history_ptr, argv, END, h_size, reps, aeps, t_adaptive_size_ptr, t_adaptive_ptr, &numbers, functionPFR3)
        printf("step size number is: %d\n", numbers)
        if numbers >= MAX_STEPS:
            print(f"Warning! Please make sure MAX_STEPS > {numbers}")
        return t_history[:numbers], X_history[:numbers * N], t_adaptive_size[:numbers], t_adaptive[:numbers]
    elif model == 6:
        fill_array(X_ptr, x0, <int>argv[9])
        compute_PFR(t_history_ptr, X_ptr, X_history_ptr, argv, END, h_size, reps, aeps, t_adaptive_size_ptr, t_adaptive_ptr, &numbers, functionPFR1)
        printf("step size number is: %d\n", numbers)
        if numbers >= MAX_STEPS:
            print(f"Warning! Please make sure MAX_STEPS > {numbers}")
        return t_history[:numbers], X_history[:numbers * N], t_adaptive_size[:numbers], t_adaptive[:numbers]
    elif model == 7:
        compute(t_history_ptr, x_history_1_ptr, x_history_2_ptr, x_history_3_ptr, x, argv, END, h_size, reps, aeps, t_adaptive_size_ptr, t_adaptive_ptr, &numbers, functionstability)
    
        
    cdef double end = time.perf_counter()
    cdef double Times = end - start
    printf("Time: %f seconds\n", Times)
    printf("step size number is: %d\n", numbers)
    if numbers >= MAX_STEPS:
        printf("Warning! Please make sure MAX_STEPS > numbers")
    cdef cnp.ndarray[double, ndim=2] x_history = np.empty((numbers, 3), dtype=np.float64)
    cdef double* x_history_ptr = <double*>x_history.data
    cdef int i = 0
    for i in range(numbers):
        x_history_ptr[3*i    ] = x_history_1[i]
        x_history_ptr[3*i + 1] = x_history_2[i]
        x_history_ptr[3*i + 2] = x_history_3[i]
    
    return t_history[:numbers], x_history, t_adaptive_size[:numbers], t_adaptive[:numbers]