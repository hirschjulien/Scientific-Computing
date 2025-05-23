function [W,Tw,dW] = StdWienerProcess(T,N,nW,Ns,seed)
% StdWienerProcess    Ns realizations of a standard Wiener process
%
% Syntax: [W,Tw,dW] = StdWienerProcess(T,N,Ns,seed)
%   W   : Standard Wiener process in [0,T]
%   Tw  : Time points
%   dW  : White noise used to generate the Wiener process
%
%   T   : Final time
%   N   : Number of intervals
%   nW  : Dimension of W(k)
%   Ns  : Number of realizations
%   seed: To set the random number generator (optional)

if nargin == 4
    rng(seed);
end

dt = T/N;
dW = sqrt(dt) * randn(nW,N,Ns);
W = [zeros(nW,1,Ns) cumsum(dW,2)];
Tw = 0:dt:T;
end
