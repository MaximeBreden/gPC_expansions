clear variables 
close all
clc

% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%% INITIALIZATION

%parameters of the model
para.sigma = 10; %not random 
para.beta = 8/3; %not random
rho_bar = 28; %we consider rho = rho_bar + delta*p with p variable
delta = 10;

PC_type = 'Legendre'; % Choice of basis for the gPC expansion
mu = 20; % Only meaningful for Gegenbauer polynomials
% The implemented expansions are:
% 'Legendre'
% 'Chebyshev' (of the first kind)
% 'Chebyshev2' (of the second kind)
% 'Gegenbauer' (of parameter mu specified above)
% 'Taylor' (corresponding to the monomial basis)
%
% /!\ The normalization is not the "standard" one for Chebyshev2 and
% Gegenbauer (see the paper for more details)

K = 100; % Truncation level in Fourier
N = 15; % Truncation level for the gPC expansion

str = ['OrbitLorenz_rho',num2str(rho_bar),'_delta',num2str(delta),'_',PC_type,'.mat'];
if exist(str,'file') == 2
    load(str,'X')
else
    error('There is no precomputed orbit for this choice of parameters and gPC basis, but you can try to get one using script_LorenzExplore.m')
end

X = resize(X,K,N);
para.Xphase = X(2:end,:);


%% REFINEMENT OF THE APPROXIMATE SOLUTION USING NEWTON'S METHOD

show = 1;
it_max = 20;
tol = 10^-10;
para.alpha_mat = reshape(Tens_for_prod(N,PC_type,mu),[N^2,N]); %linearization coefficients
para.rho = [rho_bar delta zeros(1,N-2)];
para.symmetrize = @symmetrize;

fprintf("\nRefinement of the stored solution using Newton's method\n")
X = Newton(X, @F_Lorenz, @DF_Lorenz, para, it_max, tol, show);

% %% PLOTS
% figure(10)
% clf
% plot_Lorenz(eval_PC(X,-1,PC_type,mu),10000,'-r',10,11)
% plot_Lorenz(eval_PC(X,0,PC_type,mu),10000,'-g',10,11)
% plot_Lorenz(eval_PC(X,1,PC_type,mu),10000,'-b',10,11)
% figure(10)
% legend('p=-1','p=0','p=1')

%% PREVALIDATION (without interval arithmetic)
para.PC_type = PC_type;
para.mu = mu;
para.nu = 1; % first weight for the norm (\ell^1_\nu)
para.eta = 1; % second weight for the norm (\ell^1_\eta)
[rmin, rmax, a] = proof_Lorenz(X, para);

%% RIGOROUS PROOF (with interval arithmetic)
if rmin < rmax
    if exist('intval','file')
        iX = intval(X);
        ia = intval(a);
        ipara = para;
        ipara.Xphase = intval(para.Xphase);
        ipara.nu = intval('1');
        ipara.eta = intval('1');
        ipara.mu = intval(mu);
        [irmin,irmax] = proof_Lorenz(iX, ipara, ia);
    else
        fprintf("\nYou need Intlab in order to run the rigorous proof\n")
    end
end



