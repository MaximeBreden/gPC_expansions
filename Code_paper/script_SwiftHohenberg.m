clear variables 
close all
clc 

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% INITIALIZATION
branch = 1; %Change 1 to 2 to select the second branch discussed in the paper
PC_type = 'Chebyshev';
str = ['SH_branch',num2str(branch),'_',PC_type,'.mat'];
load(str,'u','para') % an approximate branch of solutions obtained using script_SwiftHohenberg_Explore.m
[K, N] = size(u)

% %%% If you want to take a larger K or N, uncomment and modify the code below
% new_K = 20;
% new_N = 50;
% new_u = zeros(new_K,new_N);
% new_u(1:K,1:N) = u;
% u = new_u;
% K = new_K;
% para.rho = [para.rho, zeros(1,new_N-N)];
% para.alpha_mat = reshape(Tens_for_prod(new_N,para.PC_type,para.mu),[new_N^2,new_N]); %linearization coefficients

%% REFINEMENT OF THE APPROXIMATE SOLUTION USING NEWTON'S METHOD
show = 1;
it_max_Newton = 50;
tol_Newton = 10^-10;
fprintf("\nRefinement of the stored solution using Newton's method\n")
u = Newton(u, @F_SH, @DF_SH, para, it_max_Newton, tol_Newton, show);

%% PREVALIDATION (without interval arithmetic)
para.nu = 1; % first weight for the norm (\ell^1_\nu)
para.eta = 1; % second weight for the norm (\ell^1_\eta)
rstar = 0.01;
[rmin, rmax, a] = proof_SH(u, para, rstar);

%% RIGOROUS PROOF (with interval arithmetic)
if exist('intval','file')
    iu = intval(u);
    ia = intval(a);
    ipara = para;
    ipara.nu = intval('1');
    ipara.eta = intval('1');
    ipara.mu = intval(para.mu);
    [irmin,irmax] = proof_SH(iu, ipara, rstar, ia);
else
    fprintf("\nYou need Intlab in order to run the rigorous proof\n")
end

%% PLOTS OF THE BRANCH

% On the bifurcation diagrams
nb_p = 200;
val_delta = linspace(-1,1,nb_p);
rho_bar = para.rho(1);
delta = para.rho(2);
Tab_rho = rho_bar+val_delta*delta;
Tab_u0 = 0*Tab_rho;
Tab_l2 = 0*Tab_rho;
for ind = 1:nb_p
    up = eval_PC(u,val_delta(ind),para.PC_type,para.mu);
    Tab_u0(ind) = sum( [1, 2*ones(1,K-1)]' .* up );
    Tab_l2(ind) = sqrt( sum(up.^2) );
end

open('bif_diag_SH_u0.fig')
hold on
plot(Tab_rho, Tab_u0, 'b', 'Linewidth', 4)
xlabel('$\rho$', 'Interpreter', 'latex')
ylabel('$u(0)$', 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight

open('bif_diag_SH_l2.fig')
hold on
plot(Tab_rho, Tab_l2, 'b', 'Linewidth', 4)
xlabel('$\rho$', 'Interpreter', 'latex')
ylabel('$\Vert u \Vert_{\ell^2}$', 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight

% Several solutions along the branch
nb_p = 20;
val_delta = linspace(-1,1,nb_p);
%Legend = cell(nb_p,1);
nb_pts_plot = 100;

start = [0, 1, 1];
finish = [0, 0, 1];
% start = [0, 1, 0];
% finish = [0, 0.5, 0.2];
colors_p = [linspace(start(1),finish(1),nb_p)', linspace(start(2),finish(2),nb_p)', linspace(start(3),finish(3),nb_p)'];

figure
clf
hold on
lw = 1.5;
for ind = 1:nb_p
    plot_SH(eval_PC(u,val_delta(ind),para.PC_type,para.mu), para.L, nb_pts_plot, colors_p(ind,:), lw)
    %Legend{ind} = strcat('\rho = ', num2str(rho_bar + delta*val_delta(ind)));
end
%legend(Legend)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$u$', 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight




