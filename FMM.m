function [GreenKernel] = FMM(x,y,r0,L,k0)
% number of quadrature points 
nbtheta = L + 1;
nbphi = 2*L+1;
%r0 = x0 - y0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% quadrature on theta: L+1 points in [0,pi]
[w_theta,x_theta,~] = GaussLeg(nbtheta); % generate the Gauss-Legendre quadrature
% points for theta
theta_l = acos(x_theta);
points_theta_rep = repmat(theta_l,[1,nbphi]);
points_theta_rep = sort(points_theta_rep);
% weights for theta
w_theta_rep = [];
for j=1:nbtheta
    w_theta_rep  = [ w_theta_rep   ; w_theta(j)  *ones(nbphi, 1) ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% quadrature on phi: 2L+1 points in [0,2pi] 
m = 0; l = 0;
[Int_value] = Int_spherical_harm(x_theta,w_theta,l,m);  % test the integration over the spherical harmonics
% points for phi
i = 0:2*L;
phi_i = 2*pi.*i/nbphi;
points_phi_rep = repmat(phi_i,[1,nbtheta]);
% weights for phi
w_phi = 2*pi / nbphi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% transform quadrature points into cartesian coordinates
[Sx, Sy, Sz] = sph2cart(points_phi_rep', (points_theta_rep-pi/2)',1); 
PointsQuad = [Sx, Sy, Sz];
WeightsQuad = w_theta_rep.*w_phi;
% figure(4)
% myfig = plot3(Sx,Sy,Sz, '--ob') ;% plot the quadrature points
% set(myfig,'LineWidth',0.001)
% view(2);
% axis('equal');
% xlabel('$x$','Interpreter','latex','FontSize',24);
% ylabel('$y$','Interpreter','latex','FontSize',24);
% zlabel('$z$','Interpreter','latex','FontSize',24);
% grid on;
% %xt = get(gca, 'XTick');
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize', 24);
% % ajouter eventuellement un titre
% %title(titre);
% view(3)
% x0=50;
% y0=50;
% width=600;
% height=500;
% set(gcf,'units','points','position',[x0,y0,width,height])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute the plane wave expansion
% Ceps = 7.5; D = 0.3*lambda; L_advice = sqrt(3)*k0*D + Ceps*log(sqrt(3)*k0*D + pi);
[val] = G_onde_plane(PointsQuad,r0,L,k0); % leg pol remove the loop
%tot_exp = WeightsQuad.* exp(-1i*k0*PointsQuad*(y-y0)') .* val.' .* exp(1i*k0*PointsQuad*(x-x0)') ;
tot_exp = WeightsQuad.* exp(1i*k0*PointsQuad*(x-y+r0)') .* val.' ;
GreenKernel = sum(tot_exp); % do a big function for FMM
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%