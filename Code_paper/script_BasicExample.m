clear variables 
close all
clc

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% INITIALIZATION
g = [2, -1, 0, 0, 0, -0.5]'; 
syms t
figure(1)
fplot(legendreP(0:5,t)*g, [-1,1], '--r','LineWidth',2)
hold on

N = 6; %number of gPC coefficients that we use for x
PC_type = 'Legendre';
if length(g) > N
    error("N is too small")
else
    g = [g; zeros(N-length(g),1)];
end

%% COMPUTATION OF AN APPROXIMATE SOLUTION
fprintf("\nComputation of an approximate solution\n")
x = zeros(N,1);
x(1) = sqrt(g(1)); %rough initialization for Newton's method
para.g = g;
para.alpha_mat = reshape(Tens_for_prod(N,PC_type),[N^2,N]); %the linearization coefficients

show = 0; %Put to 1 (resp. 0) to show (resp. hide) the errors during Newton's iterations
it_max_Newton = 50;
tol_Newton = 10^-12;
x = Newton(x, @F_Example, @DF_Example, para, it_max_Newton, tol_Newton, show);

%% PREVALIDATION (without interval arithmetic)
para.PC_type = PC_type;
para.mu = 0; %meaningful only for Gegenbauer expansion of for using interval arithmetic
para.eta = 1.0; %the weight in the gPC "direction" (\ell^1_\eta)
[rmin, rmax, a] = proof_Example(x, para);

%% RIGOROUS PROOF (with interval arithmetic)
if exist('intval','file') && not(isnan(rmin))
    ix = intval(x);
    ia = intval(a);
    ipara = para;
    ipara.g = intval(g);
    ipara.eta = intval('1.0');
    ipara.mu = intval(0);
    [irmin,irmax] = proof_Example(ix, ipara, ia);
else
    fprintf("\nYou need Intlab in order to run the rigorous proof\n")
end

%% PLOT OF THE SOLUTION
nb_pts = 500;
fig = 1;
lw = 2;
plot_PC(x, PC_type, 0, 'b', nb_pts, fig, lw)
xlabel('$p$', 'Interpreter', 'latex')
legend('$g$', '$\bar{x}$', 'Interpreter', 'latex')
set(gca,'FontSize',15)
axis tight
drawnow

fprintf("\nThe approximate mean value of x is %14.13f\n",x(1))
