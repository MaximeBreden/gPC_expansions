clear variables 
close all
clc 

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% INITIALIZATION
load('SH_275_3.mat','u','para') % an approximate solution obtained for given parameters using the continuation script
para.rho = 3;
para.beta = 1;
K = 20; % number of Fourier modes
if K <= length(u)
    u = u(1:K);
else
    u = [u;zeros(K-length(u),1)];
end

%% REFINEMENT OF THE APPROXIMATE SOLUTION USING NEWTON'S METHOD
show = 1;
it_max_Newton = 50;
tol_Newton = 10^-10;
fprintf("\nRefinement of the stored solution using Newton's method\n")
u = Newton(u, @F_SH_deter, @DF_SH_deter, para, it_max_Newton, tol_Newton, show);

% figure
% nb_pts_plot = 100;
% plot_SH(u, para.L, nb_pts_plot, 'k')

%% COMPUTATION OF A FOURIER x gPC REPRESENTATION OF A BRANCH
% This is done with Newton's method, initialized with delta=0, and possibly using continuation in delta.  
delta = 1; % The varying parameter is of the form rho_bar + delta*p, where p varies in [-1,1]
Tab_delta = linspace(0,delta,5); % If delta is too large for Newton's method to converge, add more intermediate values in Tab_delta 

PC_type = 'Chebyshev'; % Choice of basis for the gPC expansion
% The implemented expansions are:
% 'Legendre'
% 'Chebyshev' (of the first kind)
% 'Chebyshev2' (of the second kind)
% 'Gegenbauer' (of parameter mu specified below)
% 'Taylor' (corresponding to the monomial basis)
%
% /!\ The normalization is not the "standard" one for Chebyshev2 and
% Gegenbauer (see the paper for more details)

mu = 20; % Only meaningful for Gegenbauer polynomials
N = 10; % Truncation level for the gPC expansion

u = [u, zeros(K,N-1)];
u_old = u;
para.alpha_mat = reshape(Tens_for_prod(N,PC_type,mu),[N^2,N]); %linearization coefficients
rho_bar = para.rho;
show = 1;
fprintf("\nNumerical continuation in delta:\n")
for delta = Tab_delta
    fprintf("\ndelta = %f\n",delta)
    para.rho = [rho_bar delta zeros(1,N-2)]; %assumes that N >= 2
    u = Newton(u, @F_SH, @DF_SH, para, it_max_Newton, tol_Newton, show);
end

%% PREVALIDATION (without interval arithmetic)
para.PC_type = PC_type;
para.mu = mu;
para.nu = 1; % first weight for the norm (\ell^1_\nu)
para.eta = 1; % second weight for the norm (\ell^1_\eta)
rstar = 0.01;
[rmin, rmax, a] = proof_SH(u, para, rstar);

%% RIGOROUS PROOF (with interval arithmetic)
if rmin < rmax
    if exist('intval','file')
        iu = intval(u);
        ia = intval(a);
        ipara = para;
        ipara.nu = intval('1');
        ipara.eta = intval('1');
        ipara.mu = intval(mu);
        [irmin,irmax] = proof_SH(iu, ipara, rstar, ia);
    else
        fprintf("\nYou need Intlab in order to run the rigorous proof\n")
    end
end

%% PLOTS OF THE BRANCH

% On the bifurcation diagrams
nb_p = 200;
val_delta = linspace(-1,1,nb_p);
Tab_rho = rho_bar+val_delta*delta;
Tab_u0 = 0*Tab_rho;
Tab_l2 = 0*Tab_rho;
for ind = 1:nb_p
    up = eval_PC(u,val_delta(ind),PC_type,mu);
    Tab_u0(ind) = sum( [1, 2*ones(1,K-1)]' .* up );
    Tab_l2(ind) = sqrt( sum(up.^2) );
end

open('bif_diag_SH_u0.fig')
hold on
plot(Tab_rho, Tab_u0, 'b', 'Linewidth', 4)
drawnow

open('bif_diag_SH_l2.fig')
hold on
plot(Tab_rho, Tab_l2, 'b', 'Linewidth', 4)
drawnow

% Several solutions along the branch
nb_p = 10;
val_delta = linspace(-1,1,nb_p);
Legend = cell(nb_p,1);
nb_pts_plot = 100;

start = [0, 1, 1];
finish = [0, 0, 1];
colors_p = [linspace(start(1),finish(1),nb_p)', linspace(start(2),finish(2),nb_p)', linspace(start(3),finish(3),nb_p)'];

figure
clf
hold on
lw = 1.5;
for ind = 1:nb_p
    plot_SH(eval_PC(u,val_delta(ind),PC_type,mu), para.L, nb_pts_plot, colors_p(ind,:), lw)
    Legend{ind} = strcat('\rho = ', num2str(rho_bar + delta*val_delta(ind)));
end
legend(Legend)

%% EXTENSION TO THE CASE WITH 2 PARAMETERS
N1 = N;
N2 = 10;
u2 = zeros(K,N1,N2);
u2(:,:,1) = u;
para2 = para;

beta_bar = para.beta;
delta_beta = 0.5; % The second varying parameter is of the form beta_bar + delta_beta*p2, where p2 varies in [-1,1]
Tab_delta_beta = linspace(0,delta_beta,5); % If delta_beta is too large for Newton's method to converge, add more intermediate values in Tab_delta 
PC_type2 = 'Chebyshev';
mu2 = 20; 
para2.alpha_mat = { para.alpha_mat, reshape(Tens_for_prod(N2,PC_type2,mu2),[N2^2,N2]) };

%%% Numerics to get an approximate manifold of solutions
fprintf("\nNumerical continuation in delta_beta:\n")
for delta_beta = Tab_delta_beta
    fprintf("\ndelta_beta = %f\n",delta_beta)
    para2.beta = [beta_bar delta_beta zeros(1,N2-2)]; %assumes that N2 >= 2
    u2 = Newton(u2, @F_SH_2para, @DF_SH_2para, para2, it_max_Newton, tol_Newton, show);
end


%%% Prevalidation (without interval arithmetic)
para2.PC_type = {PC_type, PC_type2};
para2.mu = [mu mu2];
para2.eta = [1 1];
rstar = 0.01;
[rmin2, rmax2, a2] = proof_SH_2para(u2, para2, rstar);

%%% Rigorous proof (with interval arithmetic)
if rmin2 < rmax2
    if exist('intval','file')
        iu2 = intval(u2);
        ia2 = intval(a2);
        ipara2 = para2;
        ipara2.nu = intval('1');
        ipara2.eta = [intval('1'),intval('1')];
        ipara2.mu = intval(para2.mu);
        [irmin2,irmax2] = proof_SH_2para(iu2, ipara2, rstar, ia2);
    else
        fprintf("\nYou need Intlab in order to run the rigorous proof\n")
    end
end

%%% Plots
fprintf("\nPlotting the manifold of solutions...\n")
nb_p1 = 20;
nb_p2 = 20;
val_p1 = linspace(-1,1,nb_p1);
val_p2 = linspace(-1,1,nb_p2);
[RHO, BETA] = meshgrid( rho_bar + delta*val_p1, beta_bar + delta_beta*val_p2 );
u_eval_PC = zeros(K,nb_p2,nb_p1);
for i1 = 1:nb_p1
    for i2 = 1:nb_p2
        p = [val_p1(i1), val_p2(i2)];
        u_eval_PC(:,i2,i1) = eval_PC_2para(u2, p, para2.PC_type, para2.mu);
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



