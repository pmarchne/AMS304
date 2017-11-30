function [val] = Int_spherical_harm(xi,wi,l,m)
% 0 <= l <= 2L(=m)
% -l < m < l

Const = (-1)^(m + abs(m))*sqrt(   (2*l+1)*factorial(l-abs(m)) / (4*pi*factorial(l+abs(m)))   ); % harmoniques spheriques normalisees

L = legendre(l,xi);
L_leg = L(abs(m+1),:);
sumLeg = wi*(L_leg)';

%%%%%%%%%%%%%%%%%%%%%%%% quadrature uniforme %%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = l-1 ;
j = 0:J;
wi_phi = 2*pi ./ l*ones(l,1);
phi_i = 2*pi.*j/l;
quad_uni = sum(wi_phi'*(exp(1i*m.*phi_i))');

if (l==0)
    quad_uni = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

val = Const * sumLeg * quad_uni;


end

