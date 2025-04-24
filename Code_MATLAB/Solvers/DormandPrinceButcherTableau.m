function [c,AT,b] = DormandPrinceButcherTableau()
% Dormant prience butcher tableau
% The Dormant prince method is of order 5, error estimation order 4,
% therefore DOPRI5(4)
c=[0,       1/5, 3/10,      4/5,     8/9,        1,     1];
b=[35/384,  0,   500/1113,  125/192, -2187/6784, 11/84, 0];
%A=zeros(length(c),length(b));
A = [0,          0,           0,          0,        0,           0,     0;
     1/5,        0,           0,          0,        0,           0,     0;
     3/40,       9/40,        0,          0,        0,           0,     0;
     44/45,      -56/40,      32/9,       0,        0,           0,     0;
     19372/6561, -25360/2187, 64448/6561, -212/729, 0,           0,     0;
     9017/3168,  -355/33,     46732/5247, 49/176,   -5103/18656, 0,     0;
     35/384,     0,           500/1113,   125/192,  -2187/6784,  11/84, 0];

AT=A';