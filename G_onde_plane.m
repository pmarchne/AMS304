function [val] = G_onde_plane(s,r0,L,k0)
%Q=(L+1)*(2*L+1);
% s is a Q*3 vector
% r0 is a 1*3 vector

Nquad = length(s);
r0_mat = repmat(r0,[Nquad,1]);

kr0 = k0*norm(r0,2);
Const = (1i * k0) / (4*pi);

norms =  sqrt(s(:,1).^2 + s(:,2).^2 + s(:,3).^2);
cos_theta = dot(s,r0_mat,2) ./ (norm(r0,2).*norms);

ind = 0:L;
hp1 = sphbessel(ind, (kr0) );

% All this is aimed to avoid generating the legendre polynomial matrix
% several times
name = 'LegPoyCosMat.mat';
if (~exist(name))
    legendrePmat = zeros(L+1,Nquad);
    for q = 0:L
        legendrePmat(q+1,:) = legendreP(q, cos_theta) ;
    end
    save(name,'legendrePmat');
else
    load(name,'legendrePmat');
    L_old = size(legendrePmat,1) -1;
    if(L_old ~= L)
        legendrePmat = zeros(L+1,Nquad);
        for q = 0:L
          legendrePmat(q+1,:) = legendreP(q, cos_theta) ;
        end 
        save(name,'legendrePmat');
        load(name,'legendrePmat');
    end
end

val = Const.*(2.*ind+1).*(1i.^ind).*hp1*legendrePmat;
end

