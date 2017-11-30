function affichemaillage(nom_maillage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% affichemaillage: 
% pour visualiser un maillage triangulaire 2D
%
% SYNOPSIS affichemaillage(nom_maillage, titre)
%          
% INPUT  * nom_maillage : racine du fichier de maillage .msh (string) 
%        * titre (optionel) un titre (string)
%
% OUTPUT une fenetre graphique
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% control on the input args
%if (nargin<2), titre = ''; end;

% lecture du fichier nom_maillage.amdba
[coor,tri,~]=read_meshfile(nom_maillage);
%visualisation du maillage
figure;
hold on

% maillage
h = trimesh(tri,coor(:,1),coor(:,2),coor(:,3),'Edgecolor','blue');
set(h,'LineWidth',0.001)
view(2);
axis('equal');
xlabel('$x$','Interpreter','latex','FontSize',24);
ylabel('$y$','Interpreter','latex','FontSize',24);
zlabel('$z$','Interpreter','latex','FontSize',24);
grid on;
%xt = get(gca, 'XTick');
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 24);
% ajouter eventuellement un titre
%title(titre);
view(3)
x0=50;
y0=50;
width=600;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])

hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
