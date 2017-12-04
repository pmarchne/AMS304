%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A routine for the implementation of a single-level Fast Multipole Method
% (FMM) over the unit sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
cd('/home/philippe/Documents/MATLAB/AMS304/Projet');
clear all; 
clc; close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import meshes
addpath('/home/philippe/Documents/MATLAB/AMS304/Projet/Mesh');
MeshNames = {'sphere_mesh_0','sphere_mesh_1','sphere_mesh_2','sphere_mesh_3','sphere_mesh_4','sphere_mesh_5'};
Mesh0Sphere = MeshNames{1,1}; k0 = 4; % associated wavenumber to get ~10 points per wavelength
Mesh1Sphere = MeshNames{1,2}; k1 = 8;
% affichemaillage(Mesh1Sphere); % uncomment to plot the mesh
Mesh2Sphere = MeshNames{1,3}; k2 = 16;
Mesh3Sphere = MeshNames{1,4}; k3 = 40;

k = [k0, k1, k2]; 
%k = [2, 4, 8 , 16 , 32, 64, 128, 256, 512, 1024, 2048];
lambda = 2*pi./k; 

% select the relevant wavenumber
k0 = k(1);
lambda0 = 2*pi./k0;

% import the nodes and triangle from teh mesh file
[NodesCoor,NodesTri,~]=read_meshfile(Mesh0Sphere);
NumberOfNodes = length(NodesCoor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classic technique complexity
GreenMat = GenGreenMat(NodesCoor,k0); % compute the Green matrix with classic technique

rho = rand(NumberOfNodes,1);
V = GreenMat*rho; % Result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% complexity test: classic matrix-vector product technique
% complexity test: uncomment to use it

% % N = [400, 800 , 1600 , 3200, 6400, 12800, 25600];
% %ElapsedTime = zeros(length(N),1);
% % count = 1;
% for p = 1:length(k)
%  %for p = N
%      [NodesCoor,~,~]=read_meshfile(MeshNames{1,p});
%      N(p) =  length(NodesCoor);
%      %tic;
% %      rho = rand(p,1);
% %      TmpMat = rand(p,p);
%      GreenMatTest = GenGreenMat(NodesCoor,k(p));
%      tic;
%      %ElapsedTime(p) = toc;
%      ResultTest = GreenMat*rho;
%      %ResultTest = TmpMat*rho;
%      ElapsedTime(p) = toc;
%      %count = count + 1;
%  end
% figure(2)
% loglog(N./k,ElapsedTime,'--+r','linewidth',1.5) % plot complexity graph
% hold on;
% loglog(N./k,8e-8*(N./k).^2,'-ok','linewidth',1.5) 
% xlabel('$N$ ', 'Interpreter' , 'latex','FontSize', 24);
% ylabel('time','Interpreter','latex','FontSize',24);
% %xlim([3 20]);
% %ylim([0.01 20]);
% grid on;
% legend({'Classic matrix-vector product','$\mathcal{O}(N^2)$'},'Interpreter','latex')
% 
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize', 20);
% 
% x0=50;
% y0=50;
% width=600;
% height=500;
% set(gcf,'units','points','position',[x0,y0,width,height])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fast multipole method

% test conditions
%x0=[0, 0, 0];
%x=x0-0.1*ones(1, 3);
y = [0,5,1]; % integration point
y0= [98,0.2,1];
%y=y0+0.03*ones(1, 3);
x = [-3,5,40]; % collocation point
x0 = [3,5,2.5];
k0 =10; lambda = 2*pi/k0;
r0 = x0 - y0;
r = x - y - r0;
% condition: abs(r0) > abs(r)
tmp = abs(r-r0)/abs(r0) % < 1, 0.8944
testr0 = norm(r0,2)/(0.6*lambda) % doit etre grand par rapport a 1
if(testr0 <= 1)
    disp('warning : r0 norm should be higher than 0.6 lambda !');
end

if(abs(r0) <= abs(r))
    disp('warning : r0 should be higher than r in magnitude !');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% complexity for L points Gauss-Legendre quadrature
% uncomment to use it
% 
% p = 1;
% N = [80,100,200,300,400,500,1000,2000];
% for n = N
%     [wi,xi,time(p)] = GaussLeg(n);
%     p = p+1;
% end
% figure(3)
% loglog(N,time,'--+r','linewidth',1.5) % plot complexity graph
% hold on;
% loglog(N,0.00000001*N.^2.3,'-k','linewidth',1) 
% xlabel('$N$ ', 'Interpreter' , 'latex','FontSize', 24);
% ylabel('time','Interpreter','latex','FontSize',24);
% % xlim([1 1000]);
% % ylim([0.005 1]);
% legend({'Gauss-Leg','$\mathcal{O}(N^{2.3})$'},'Interpreter','latex')
% grid on;
% 
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize', 20);
% 
% x0=50;
% y0=50;
% width=600;
% height=500;
% set(gcf,'units','points','position',[x0,y0,width,height])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Implementation of the plane wave expansion
%L = 10; % number of quadrature points 
%[rep] = FMM(x,y,x0,y0,L,k0);
%% Exact solution
 
reponse_ex = exp(1i * k0 * norm(x-y,2))./ (4*pi*norm(x-y,2));
%err = rep - reponse_ex;
% plot err function of L for 3 configurations of (x0,y0) and (x,y)
counter = 1; 
L = 10:10:90;
for l = L
    rep(counter) = FMM(x,y,x0,y0,l,k0);
    err(counter) = rep(counter) - reponse_ex;
    counter = counter + 1;
end
figure(4)
loglog(L,abs(err),'--+r','linewidth',1.5) % plot complexity graph 
xlabel('Truncation parameter $L$ ', 'Interpreter' , 'latex','FontSize', 24);
ylabel('FMM error','Interpreter','latex','FontSize',24);
% xlim([1 1000]);
% ylim([0.005 1]);
%legend({'Gauss-Leg','$\mathcal{O}(N^{2.3})$'},'Interpreter','latex')
grid on;

set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize', 20);

x0=50;
y0=50;
width=600;
height=500;
set(gcf,'units','points','position',[x0,y0,width,height])


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Octree decomposition
% 
% maxsize = 0.3*pi/2; % cell size needs to be greatr than 0.3*lambda0
% OT = OcTree(NodesCoor,'binCapacity',length(NodesCoor),'maxSize', maxsize);        
% figure
% boxH = OT.plot;
% cols = lines(OT.BinCount);
% doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
% for i = 1:OT.BinCount
%     set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
%     doplot3(NodesCoor(OT.PointBins==i,:),'.','Color',cols(i,:))
% end
% axis image, view(3)
% 
% % interaction proche si |x0-y0|< 0.6*lambda
% pas = 0.125;
% centrex = -1+pas:pas:1-pas;
% % if inter_proch do green classic
% % if not determine corresponding x0 and y0 and do FMM