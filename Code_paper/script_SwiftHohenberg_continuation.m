clear variables 
close all
clc

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% Initialization
 
para.L = 2*pi;% the domain is (0,L)
para.beta = 1;
para.rho = 0.3;
K = 50; % number of Fourier modes
u = zeros(K,1);

%% Numerical continuation

pt0 = [para.rho; u];
vect0 = [1; zeros(K,1)];
ds = 1e-3; % step-size for the pseudo arclength continuation
depth = 5; % maximal depth of the tree (i.e. maximal number of bifurcations from the starting point)
nb_pts_max = 10000; % maximal number of points on a branch (between 2 successive bifurcations)
rho_max = 4.5; % maximal value of rho
it_max_Newton = 50;
tol_Newton = 10^-10;
[tab_sol, tab_bif] = continuation_SH_pal(@F_SH_deter, @DF_SH_deter_cont, pt0, para, vect0, ds, depth, nb_pts_max, rho_max, tol_Newton, it_max_Newton);  

%% Bifurcation diagrams

tab_rho = tab_sol(1,:);
tab_rho_bif = tab_bif(1,:);

tab_u0 = sum( diag([1, 2*ones(1,K-1)]) * tab_sol(2:end,:), 1);
tab_u0_bif = sum( diag([1, 2*ones(1,K-1)]) * tab_bif(2:end,:), 1);
figure(10)
clf
plot(tab_rho, tab_u0, '.k')
hold on
plot(tab_rho_bif, tab_u0_bif, '*r', 'Linewidth', 5)

tab_normu = sqrt( sum( tab_sol(2:end,:).^2, 1) );
tab_normu_bif = sqrt( sum( tab_bif(2:end,:).^2, 1) );
figure
plot(tab_rho, tab_normu, '.k')
hold on
plot(tab_rho_bif, tab_normu_bif, '*r', 'Linewidth', 5)


% %% To find a solution for a specific value of rho
% 
% wanted_rho = 2.75;
% margin = ds/2;
% indices = find( abs(tab_rho - wanted_rho) < margin );
% for j = 1:length(indices)
%     figure
%     plot_SH(tab_sol(2:end,indices(j)), para.L, 100)
%     drawnow
%     pause
% end









