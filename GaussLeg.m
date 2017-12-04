% Quadrature de Gauss-Legendre à L points. Cette quadrature est exacte
% pour les polynomes de degrés au plus 2L-1. Les poids sont déterminés par
% les racines des polynomes de Legendre.

function [wi,xi,time] = GaussLeg(L)
tic;
xi =zeros(1,L);
wi = zeros(1,L);
% ------------------------------
d = [-1,1];
I = 1:(L-1);
T = I ./sqrt(4.*I.^2 -1); T1 = [T,0] ;T2 = [0,T];
A = spdiags([T1',T2'],d,L,L);
A = full(A);
% ------------------------------
[V,D] = eig(A);
xi = (diag(D))';
wi = 2*(V(1,:).^2);
time = toc;
end

%somme des poids = 2