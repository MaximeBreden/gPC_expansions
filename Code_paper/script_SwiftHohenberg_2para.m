clear variables 
close all
clc 

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% INITIALIZATION
load('SH_manifold','u','para') % an approximate manifold of solutions obtained using script_SwiftHohenberg_Explore.m
[K, N1, N2] = size(u)

%% REFINEMENT OF THE APPROXIMATE SOLUTION USING NEWTON'S METHOD
show = 1;
it_max_Newton = 50;
tol_Newton = 10^-10;
fprintf("\nRefinement of the stored solution using Newton's method\n")
u = Newton(u, @F_SH_2para, @DF_SH_2para, para, it_max_Newton, tol_Newton, show);

%% PREVALIDATION (without interval arithmetic)
para.nu = 1; % first weight for the norm (\ell^1_\nu)
para.eta = [1,1]; % second weights for the norm (\ell^1_\eta)
rstar = 0.01;
[rmin, rmax, a] = proof_SH_2para(u, para, rstar);

%% RIGOROUS PROOF (with interval arithmetic)
if exist('intval','file')
    iu = intval(u);
    ia = intval(a);
    ipara = para;
    ipara.nu = intval('1');
    ipara.eta = [intval('1'),intval('1')];
    ipara.mu = intval(para.mu);
    [irmin,irmax] = proof_SH_2para(iu, ipara, rstar, ia);
else
    fprintf("\nYou need Intlab in order to run the rigorous proof\n")
end

%% PLOTS OF THE BRANCH

nb_p1 = 20;
nb_p2 = 20;
val_p1 = linspace(-1,1,nb_p1);
val_p2 = linspace(-1,1,nb_p2);
[RHO, BETA] = meshgrid( para.rho(1) + para.rho(2)*val_p1, para.beta(1) + para.beta(2)*val_p2 );
u_eval_PC = zeros(K,nb_p2,nb_p1);
for i1 = 1:nb_p1
    for i2 = 1:nb_p2
        p = [val_p1(i1), val_p2(i2)];
        u_eval_PC(:,i2,i1) = eval_PC_2para(u, p, para.PC_type, para.mu);
    end
end

eval_Four = ones(K,nb_p2,nb_p1);
eval_Four(2:end,:,:) = 2;
manifold_u0 = squeeze(sum(u_eval_PC .* eval_Four, 1));
manifold_l2 = squeeze(sqrt(sum(u_eval_PC .^ 2, 1)));

figure
surf(RHO,BETA,manifold_u0)
shading interp
xlabel('$\rho$', 'Interpreter', 'latex')
ylabel('$\beta$', 'Interpreter', 'latex')
zlabel('$u(0)$', 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight
drawnow

figure
surf(RHO,BETA,manifold_l2)
shading interp
xlabel('$\rho$', 'Interpreter', 'latex')
ylabel('$\beta$', 'Interpreter', 'latex')
zlabel('$\Vert u \Vert_{\ell^2}$', 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight
drawnow



