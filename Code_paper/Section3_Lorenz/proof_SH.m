function [rmin, rmax, a] = proof_SH(u, para, rstar, a)

%The last argument is optional, and will be computed inside this function
%if not given as input.

if exist('intval','file') && isintval(u(1))
    fprintf('\nRigorous validation with interval arithmetic ...\n')
else
    fprintf('\nPrevalidation without interval arithmetic ...\n')
end

%% Initialization

L = para.L;
beta = para.beta;
rho = para.rho;
nu = para.nu;
eta = para.eta;
PC_type = para.PC_type;
mu = para.mu;

[K,N] = size(u);

%% Building A
if nargin < 4
    Iconvo = zeros(K*N, K);
    for k = 1:K
        Iconvo((k-1)*N+1,k) = 1;
    end
    a = DF_SH(u, para) \ Iconvo; 
    a = permute( reshape(a, [N, K, K]), [2,3,1] ); %Tensor containing the elements a_{k,l} = a(k,l,:) representing A
end
a_K4N = cat(3, a, zeros(K, K, 3*(N-1)));

%% Preparations (zero padding to remove truncation errors, weights, linearization coefficients ...)
K2 = 2*K-1;
K3 = 3*K-2;
K5 = 5*K-4;
N2 = 2*N-1;
N3 = 3*N-2;
N4 = 4*N-3;
u_3K3N = [[u;zeros(2*K-2,N)], zeros(3*K-2,2*N-2)];
u_5K3N = [u_3K3N;zeros(2*K-2,N3)];
rho_3N = [rho zeros(1,2*N-2)];

weights_K = [1 2*ones(1,K-1)] .* nu.^(0:K-1);
weights_2K = [1 2*ones(1,K2-1)] .* nu.^(0:K2-1);
weights_3K = [1 2*ones(1,K3-1)] .* nu.^(0:K3-1);
weights_5K = [1 2*ones(1,K5-1)] .* nu.^(0:K5-1);
weights_N = eta.^(0:N-1);
weights_2N = eta.^(0:N2-1);
weights_3N = eta.^(0:N3-1);
weights_4N = eta.^(0:N4-1);
weights_KN = (tens2vect(weights_K' * weights_N))';
weights_2K2N = (tens2vect(weights_2K' * weights_2N))';
weights_3K4N = (tens2vect(weights_3K' * weights_4N))';

I_tens = zeros(K5, K5, 1);
I_tens(:) = eye(K5);
I_tens = cat(3, I_tens, zeros(K5, K5, N3-1));

alpha_tens_4N = Tens_for_prod(N4,PC_type,mu);
alpha_mat_4N = reshape(alpha_tens_4N,[N4^2,N4]);
alpha_mat_3N = reshape(alpha_tens_4N(1:N3,1:N3,1:N3),[N3^2,N3]);
alpha_mat_2N = reshape(alpha_tens_4N(1:N2,1:N2,1:N2),[N2^2,N2]);

para_3N = para;
para_3N.rho = rho_3N;
para_3N.alpha_mat = alpha_mat_3N;


Lambda_K5 = (1-((0:K5-1)*pi/L).^2).^2;
if (K*pi/L)^2 >= 1
    min_lambda_tail = (1-(K*pi/L)^2)^2;
else
    warning('K may not be large enough')
    K_lim = ceil( L/pi );
    min_lambda_tail = min( (1-((K:K_lim)*pi/L).^2).^2 );
end
%% Y & Z1_finite
%AF(u)
F_3K3N_tens = F_SH(u_3K3N, para_3N);
F_K4N_vect = tens2vect( [F_3K3N_tens(1:K,:), zeros(K,N-1)] ); % finite part which will be multiplied by A later on
AF_3K4N_tens = zeros(K3,N4);
if exist('intval','file') && isintval(u(1))
    AF_3K4N_tens = intval(AF_3K4N_tens);
end  
AF_3K4N_tens(K+1:end,:) = diag(1./Lambda_K5(K+1:K3)) * [F_3K3N_tens(K+1:end,:), zeros(2*K-2,N-1)]; % tail part

%ADF(u)
DF_53K3N_tens = DF_SH_tens(u_5K3N, para_3N); %% /!\ We only need the columns up to k = 3K-3, but these columns have non zero entries up to the row 5K-5
DF_53K3N_tens = DF_53K3N_tens(:,1:K3,:);
DF_13K3N_mat = reshape( permute( DF_53K3N_tens(1:K,:,:), [3,1,2] ), [N3*K,K3]); % finite part of DF which will be multiplied by A later on
ADF_53K3N_tens = zeros(K5,K3,N3);
if exist('intval','file') && isintval(u(1))
    ADF_53K3N_tens = intval(ADF_53K3N_tens);
end  
ADF_53K3N_tens(K+1:end,:,:) = reshape( -diag(1./Lambda_K5(K+1:end)) * reshape( DF_53K3N_tens(K+1:end,:,:), [K5-K,N3*K3] ), [K5-K,K3,N3] ); % tail part

for k = 1:K % multiplication by a for the finite part, row by row (in k) for A    
    %AF(u)
    A_k4N = convo_PC_mat( squeeze(a_K4N(k,:,:)), alpha_mat_4N, 'column');
    AF_3K4N_tens(k,:) = A_k4N * F_K4N_vect;
    %ADF(u)
    A_k3N = convo_PC_mat( squeeze(a_K4N(k,:,1:N3)), alpha_mat_3N, 'column');
    ADF_53K3N_tens(k,:,:) = permute( reshape( A_k3N * DF_13K3N_mat, [N3,1,K3] ), [2,3,1] );
end

Y = weights_3K4N * abs(tens2vect(AF_3K4N_tens));
fprintf('\nY = %g\n',i2f(Y))

B_tens = I_tens(:,1:K3,:) - ADF_53K3N_tens;
Z1_finite = norm_op_tens(B_tens, weights_3N', weights_5K, weights_3K);
fprintf('\nZ1_finite = %g\n',i2f(Z1_finite))

%% Z1_tail
u_2K2N = [[u;zeros(K-1,N)], zeros(2*K-1,N-1)];
uu_2K2N = convo_cos_PC(u_2K2N,u_2K2N,alpha_mat_2N);
rho_2K2N = [[rho, zeros(1,N-1)]; zeros(2*K-2,2*N-1)];

Z1_tail = weights_2K2N * abs(tens2vect(rho_2K2N-3*beta*uu_2K2N)) / min_lambda_tail;
fprintf('\nZ1_tail = %g\n',i2f(Z1_tail))
            
Z1 = max(Z1_finite,Z1_tail);
            
%% Z2
norm_A_finite = norm_op_tens(a, weights_N', weights_K);

Z2 = max(norm_A_finite, 1/min_lambda_tail ) * 6*abs(beta) * ( weights_KN*abs(tens2vect(u)) + rstar );
fprintf('\nZ2 = %g\n',i2f(Z2))

%% Cheking that we have a contraction
[rmin, rmax] = iscontraction(Y, Z1, Z2, rstar);

