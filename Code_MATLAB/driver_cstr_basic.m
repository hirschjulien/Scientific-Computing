% driver_cstr_basic
% Script that runs the cstr_basic.m model
% Source: https://youtu.be/RqbnfjOF5tU?si=jEayseGCZRHYAgQi
% 

clc
clear

xa_in = .5;     % Boundary cond
xa = 1.0;       % Initial condition
t0 = 0;
tf = 5;

options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-6);
[T,X]= ode15s(@cstr_basic, [t0, tf], xa, options, xa_in);

plot(T,X)