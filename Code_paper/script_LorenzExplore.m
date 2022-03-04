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
rho_bar = 28; 
para.rho = rho_bar; %we consider rho_bar + delta*p with p variable(delta = 0 for now)

K = 100; %truncation level in Fourier

%parameters for Newton's method
show = 0;
it_max = 20;
tol = 10^-10;

%Loading a deterministic periodic orbit (there are only orbits stored for
%sigma = 10, beta = 8/3 and rho = 28 or rho = 50, but one can use 
%continuation from these to obtain deterministic periodic orbits at other 
%parameter values)
load (['AABAB_',num2str(para.rho)])
X = resize_deter(u,K);
para.Xphase = X(2:end);
% plot_Lorenz(X,1000,'r')

para.symmetrize = @symmetrize;
X = Newton(X, @F_Lorenz_deter, @DF_Lorenz_deter, para, it_max, tol,show);
%Updating the phase condition
para.Xphase = X(2:end);
% plot_Lorenz(X,1000,'b')

%% VALIDATION OF THE STARTING DETERMINISTIC ORBIT
% If this is not successful, the validation of a random orbit (which will
% include this one as a specific realisation), is already doomed to fail

para.nu = 1; % weight for the norm (\ell^1_\nu)
fprintf("\nPre-validation of the starting determinitic orbit\n")
[rmin_deter, rmax_deter] = proof_Lorenz_deter(X, para);

if not(rmin_deter<rmax_deter)
    warning('The validation of the deterministic orbit failed')
    pause
end

%% COMPUTATION OF A RANDOM ORBIT

delta = 10; % The random parameter is of the form rho_bar + delta*p, where p is random
Tab_delta = linspace(0,delta,3); % If delta is too large for Newton's method to converge, add more intermediate values in Tab_delta 

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

X = resize(X,K,N);
para.alpha_mat = reshape(Tens_for_prod(N,PC_type,mu),[N^2,N]); %linearization coefficients
show = 1;
fprintf("\n\nComputation of the gPC expansion of the random perdiodic orbit ...\n")
for delta = Tab_delta % continuation in delta, from 0 up to the desired value
    if show
        fprintf("\ndelta = %f :\n",delta)
    end
        para.rho = [rho_bar delta zeros(1,N-2)];
        fprintf("\nRunning Newton's method\n")
        X = Newton(X, @F_Lorenz, @DF_Lorenz, para, it_max, tol, show);
        %Updating the phase condition
        para.Xphase = X(2:end,:);
end

figure(10)
clf
plot_Lorenz(eval_PC(X,-1,PC_type,mu),10000,'-r',10,11)
plot_Lorenz(eval_PC(X,0,PC_type,mu),10000,'-g',10,11)
plot_Lorenz(eval_PC(X,1,PC_type,mu),10000,'-b',10,11)
figure(10)
legend('p=-1','p=0','p=1')


% str = ['OrbitLorenz_rho',num2str(rho_bar),'_delta',num2str(delta),'_',PC_type,'.mat'];
% save(str,'X')


%% PREVALIDATION (without interval arithmetic)
para.PC_type = PC_type;
para.mu = mu;
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
        ipara.mu = intval(0);
        [irmin,irmax] = proof_Lorenz(iX, ipara, ia);
    else
        fprintf("\nYou need Intlab in order to run the rigorous proof\n")
    end
end

