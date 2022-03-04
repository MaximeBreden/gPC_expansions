function [rmin, rmax, a] = proof_SH_2para(u, para, rstar, a)

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

[K,N1,N2] = size(i2f(u));
N = N1*N2;
incr_N2 = double(N2>1);

%% Building A
if nargin < 4
    Iconvo = zeros(K*N, K);
    for k = 1:K
        Iconvo((k-1)*N+1,k) = 1;
    end
    a = DF_SH_2para(u, para) \ Iconvo; 
    a = permute( reshape(a, [N2, N1, K, K]), [3,4,2,1] ); %Tensor containing the elements a_{k,l} = a(k,l,:,:) representing A
end
a_K4N = zeros(K, K, 4*N1-3, 4*N2-3+incr_N2);
if exist('intval','file') && isintval(u(1))
    a_K4N = intval(a_K4N);
end 
a_K4N(:,:,1:N1,1:N2) = a;

%% Preparations (zero padding to remove truncation errors, weights, linearization coefficients ...)

%The zero-padding used to eliminate truncation errors assumes that rho and
%beta are both 1D polynomials of degree at most 1, i.e. 
%rho = [x, x, 0, ..., 0], and similarly for beta.

K2 = 2*K-1;
K3 = 3*K-2;
K5 = 5*K-4;
N1_2 = 2*N1-1;
N1_3 = 3*N1-2;
N1_4 = 4*N1-3;
N2_2 = 2*N2-1+incr_N2;
N2_3 = 3*N2-2+incr_N2;
N2_4 = 4*N2-3+incr_N2;
N3 = N1_3*N2_3;
u_3K4N = zeros(K3,N1_4,N2_4);
if exist('intval','file') && isintval(u(1))
    u_3K4N = intval(u_3K4N);
end 
u_3K4N(1:K,1:N1,1:N2) = u;
u_5K3N = zeros(K5,N1_3,N2_3);
if exist('intval','file') && isintval(u(1))
    u_5K3N = intval(u_5K3N);
end 
u_5K3N(1:K,1:N1,1:N2) = u;
rho_3N1 = [rho zeros(1,2*(N1-1))]; 
rho_4N1 = [rho zeros(1,3*(N1-1))]; 
beta_3N2 = [beta zeros(1,2*(N2-1)+incr_N2)]; 
beta_4N2 = [beta zeros(1,3*(N2-1)+incr_N2)]; 

weights_K = [1 2*ones(1,K-1)] .* nu.^(0:K-1);
weights_2K = [1 2*ones(1,K2-1)] .* nu.^(0:K2-1);
weights_3K = [1 2*ones(1,K3-1)] .* nu.^(0:K3-1);
weights_5K = [1 2*ones(1,K5-1)] .* nu.^(0:K5-1);
weights_N1 = eta(1).^(0:N1-1);
weights_2N1 = eta(1).^(0:N1_2-1);
weights_3N1 = eta(1).^(0:N1_3-1);
weights_4N1 = eta(1).^(0:N1_4-1);
weights_N2 = eta(2).^(0:N2-1);
weights_2N2 = eta(2).^(0:N2_2-1);
weights_3N2 = eta(2).^(0:N2_3-1);
weights_4N2 = eta(2).^(0:N2_4-1);
weights_N = weights_N1' * weights_N2;
weights_2N = weights_2N1' * weights_2N2;
weights_3N = weights_3N1' * weights_3N2;
weights_4N = weights_4N1' * weights_4N2;
weights_KN = transpose( tens2vect( superkron(weights_K', reshape(weights_N,[1,N1,N2])) ) );
weights_2K2N = transpose( tens2vect( superkron(weights_2K', reshape(weights_2N,[1,N1_2,N2_2])) ) );
weights_3K4N = transpose( tens2vect( superkron(weights_3K', reshape(weights_4N,[1,N1_4,N2_4])) ) );

I_tens = zeros(K5, K5, N1_3, N2_3);
I_tens(:,:,1,1) = eye(K5);

alpha1_tens_4N1 = Tens_for_prod(N1_4,PC_type{1},mu(1));
alpha1_mat_4N1 = reshape(alpha1_tens_4N1,[N1_4^2,N1_4]);
alpha1_mat_3N1 = reshape(alpha1_tens_4N1(1:N1_3,1:N1_3,1:N1_3),[N1_3^2,N1_3]);
alpha1_mat_2N1 = reshape(alpha1_tens_4N1(1:N1_2,1:N1_2,1:N1_2),[N1_2^2,N1_2]);

alpha2_tens_4N2 = Tens_for_prod(N2_4,PC_type{2},mu(2));
alpha2_mat_4N2 = reshape(alpha2_tens_4N2,[N2_4^2,N2_4]);
alpha2_mat_3N2 = reshape(alpha2_tens_4N2(1:N2_3,1:N2_3,1:N2_3),[N2_3^2,N2_3]);
alpha2_mat_2N2 = reshape(alpha2_tens_4N2(1:N2_2,1:N2_2,1:N2_2),[N2_2^2,N2_2]);

alpha_mat_4N = {alpha1_mat_4N1,alpha2_mat_4N2};
alpha_mat_3N = {alpha1_mat_3N1,alpha2_mat_3N2};
alpha_mat_2N = {alpha1_mat_2N1,alpha2_mat_2N2};

para_3N = para;
para_3N.rho = rho_3N1;
para_3N.beta = beta_3N2;
para_3N.alpha_mat = {alpha1_mat_3N1,alpha2_mat_3N2};

para_4N = para;
para_4N.rho = rho_4N1;
para_4N.beta = beta_4N2;
para_4N.alpha_mat = {alpha1_mat_4N1,alpha2_mat_4N2};


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
F_3K4N_tens = F_SH_2para(u_3K4N, para_4N);
F_K4N_vect = tens2vect(F_3K4N_tens(1:K,:,:)); % finite part which will be multiplied by A later on
AF_3K4N_tens = zeros(K3,N1_4,N2_4);
if exist('intval','file') && isintval(u(1))
    AF_3K4N_tens = intval(AF_3K4N_tens);
end  
AF_3K4N_tens(K+1:end,:,:) = superkron(-1./Lambda_K5(K+1:K3)',ones(1,N1_4,1),ones(1,1,N2_4)) .* F_3K4N_tens(K+1:end,:,:); % tail part

%ADF(u)
DF_53K3N_tens = DF_SH_2para_tens(u_5K3N, para_3N); %% /!\ We only need the columns up to k = 3K-3, but these columns have non zero entries up to the row 5K-5
DF_53K3N_tens = DF_53K3N_tens(:,1:K3,:,:);
DF_13K3N_mat = reshape( permute( DF_53K3N_tens(1:K,:,:,:), [4,3,1,2] ), [N3*K,K3]); % finite part of DF which will be multiplied by A later on
ADF_53K3N_tens = zeros(K5,K3,N1_3,N2_3);
if exist('intval','file') && isintval(u(1))
    ADF_53K3N_tens = intval(ADF_53K3N_tens);
end  
ADF_53K3N_tens(K+1:end,:,:,:) = superkron(-1./Lambda_K5(K+1:end)',ones(1,K3,1,1),ones(1,1,N1_3,1),ones(1,1,1,N2_3)) .* DF_53K3N_tens(K+1:end,:,:,:); % tail part

extra_dim = 1;
for k = 1:K % multiplication by a for the finite part, row by row (in k) for A    
    %AF(u)
    A_k4N = convo_PC_mat( reshape(a_K4N(k,:,:,:),[K,N1_4,N2_4]), alpha_mat_4N, extra_dim, 'column');
    AF_3K4N_tens(k,:,:) = permute( reshape(A_k4N * F_K4N_vect, [N2_4,N1_4]), [2,1] );
    %ADF(u)
    A_k3N = convo_PC_mat( reshape(a_K4N(k,:,1:N1_3,1:N2_3),[K,N1_3,N2_3]), alpha_mat_3N, extra_dim, 'column');
    ADF_53K3N_tens(k,:,:,:) = permute( reshape( A_k3N * DF_13K3N_mat, [N2_3,N1_3,1,K3] ), [3,4,2,1] );
end

Y = weights_3K4N * abs(tens2vect(AF_3K4N_tens));
fprintf('\nY = %g\n',i2f(Y))

B_tens = I_tens(:,1:K3,:,:) - ADF_53K3N_tens;
Z1_finite = norm_op_tens(B_tens, weights_3N, weights_5K, weights_3K);
fprintf('\nZ1_finite = %g\n',i2f(Z1_finite))

%% Z1_tail
u_2K2N = zeros(K2,N1_2,N2_2);
if exist('intval','file') && isintval(u(1))
    u_2K2N = intval(u_2K2N);
end  
u_2K2N(1:K,1:N1,1:N2) = u;
uu_2K2N = convo_cos_PC(u_2K2N,u_2K2N,alpha_mat_2N);

rho_2K2N = zeros(K2,N1_2,N2_2);
if exist('intval','file') && isintval(u(1))
    rho_2K2N = intval(rho_2K2N);
end  
rho_2K2N(1,1:N1,1) = rho;

beta_2K2N = zeros(K2,N1_2,N2_2);
if exist('intval','file') && isintval(u(1))
    beta_2K2N = intval(beta_2K2N);
end  
beta_2K2N(1,1,1:N2) = beta;

Z1_tail = weights_2K2N * abs(tens2vect(rho_2K2N-3*convo_cos_PC(beta_2K2N,uu_2K2N,alpha_mat_2N))) / min_lambda_tail;
fprintf('\nZ1_tail = %g\n',i2f(Z1_tail))
            
Z1 = max(Z1_finite,Z1_tail);
            
%% Z2
norm_A_finite = norm_op_tens(a, weights_N, weights_K);

Z2 = max(norm_A_finite, 1/min_lambda_tail ) * 6 * weights_N2*transpose(abs(beta)) * ( weights_KN*abs(tens2vect(u)) + rstar );
fprintf('\nZ2 = %g\n',i2f(Z2))

%% Cheking that we have a contraction
[rmin, rmax] = iscontraction(Y, Z1, Z2, rstar);

