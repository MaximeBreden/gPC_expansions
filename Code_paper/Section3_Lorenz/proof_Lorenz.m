function [rmin, rmax, a] = proof_Lorenz(X, para, a)

%The last argument is optional, and will be computed inside this function
%if not given as input.

if exist('intval','file') && isintval(X(1))
    fprintf('\nRigorous validation with interval arithmetic ...\n')
else
    fprintf('\nPrevalidation without interval arithmetic ...\n')
end

%% Initialization
sigma = para.sigma;
rho = para.rho;
beta = para.beta;
nu = para.nu;
eta = para.eta;
PC_type = para.PC_type;
mu = para.mu;
alpha_mat_N = para.alpha_mat;

[K,N] = size(X);
K = (K+2)/6;

Omega = X(1,:);
x = X(2:2*K,:);
y = X(2*K+1:4*K-1,:);
z = X(4*K:6*K-2,:);

e = [1; zeros(N-1,1)];
Upsilon = convo_PC_mat( Omega, alpha_mat_N ) \ e;


%% Building A
if nargin < 3
    Iconvo = zeros((6*K-2)*N, 6*K-2);
    for k = 1:6*K-2
        Iconvo((k-1)*N+1,k) = 1;
    end
    a = DF_Lorenz(X, para) \ Iconvo; 
    a = permute( reshape(a, [N, 6*K-2, 6*K-2]), [2,3,1] ); %Tensor containing the elements a_{k,l} = a(k,l,:) representing A
end
a_K3N = cat(3, a, zeros(6*K-2, 6*K-2, 2*(N-1))); 

%% Preparations (zero padding to remove truncation errors, weights, linearization coefficients, indices sets ...)
K2 = 2*K-1;
K3 = 3*K-2;
N2 = 2*N-1;
N3 = 3*N-2;

X_2K2N = resize(X, K2, N2);
X_3K2N = resize(X, K3, N2);
rho_2N = [rho zeros(1,N-1)];

weights_K = [1 nu.^repmat(abs(-K+1:K-1), [1,3])];
weights_2K = [1 nu.^repmat(abs(-K2+1:K2-1), [1,3])];
weights_3K = [1 nu.^repmat(abs(-K3+1:K3-1), [1,3])];
weights_N = eta.^(0:N-1);
weights_2N = eta.^(0:N2-1);
weights_3N = eta.^(0:N3-1);
weights_2K3N = (tens2vect(weights_2K' * weights_3N))';
weights_singleK = nu.^abs(-K+1:K-1);
weights_singleK2N = (tens2vect(weights_singleK' * weights_2N))';

I_tens = zeros(6*K3-2, 6*K3-2, 1);
I_tens(:) = eye(6*K3-2);
I_tens = cat(3, I_tens, zeros(6*K3-2, 6*K3-2, N2-1));

alpha_tens_3N = Tens_for_prod(3*N-2,PC_type,mu);
alpha_mat_3N = reshape(alpha_tens_3N,[N3^2,N3]);
alpha_mat_2N = reshape(alpha_tens_3N(1:N2,1:N2,1:N2),[N2^2,N2]);

para_2N = para;
para_2N.rho = rho_2N;
para_2N.alpha_mat = alpha_mat_2N;

MUpsilon_2N = convo_PC_mat( [Upsilon; zeros(N-1,1)], alpha_mat_2N );
MUpsilon_3N = convo_PC_mat( [Upsilon; zeros(2*(N-1),1)], alpha_mat_3N );

% Defining indices set to easily select the appropriate modes
Ind_2K = 1:6*K2-2;
modes_2K = [0 repmat(-K2+1:K2-1, [1,3])];
Ind_Kin2K = Ind_2K(abs(modes_2K)<K); %Indices for the finite part, among indices up to 2K-2
Ind_notKin2K = Ind_2K(abs(modes_2K)>=K); %Indices for the tail part, among indices up to 2K-2
l_notKin2K = length(Ind_notKin2K);

Ind_3K = 1:6*K3-2;
modes_3K = [0 repmat(-K3+1:K3-1, [1,3])];
Ind_Kin3K = Ind_3K(abs(modes_3K)<K); %Indices for the finite part, among indices up to 3K-3
Ind_notKin3K = Ind_3K(abs(modes_3K)>=K); %Indices for the tail part, among indices up to 3K-3
Ind_2Kin3K = Ind_3K(abs(modes_3K)<K2); %Indices up to 2K-2, among indices up to 3K-3
l_notKin3K = length(Ind_notKin3K);

%% Y & Z1_finite
AntiDer_2K_tail = spdiags( 1./modes_2K(Ind_notKin2K)', 0, l_notKin2K, l_notKin2K); %The tail part of A (up to 2K-2)
AntiDer_3K_tail = spdiags( 1./modes_3K(Ind_notKin3K)', 0, l_notKin3K, l_notKin3K); %The tail part of A (up to 3K-3)

%AF(X)
F_2K2N_mat = F_Lorenz(X_2K2N, para_2N);
F_K3N_vect = tens2vect( [F_2K2N_mat(Ind_Kin2K,:), zeros(6*K-2,N-1)] ); % finite part of F which will be "multiplied" by A later on
AF_2K3N_mat = zeros(6*K2-2,N3);
if exist('intval','file') && isintval(X(1))
    AF_2K3N_mat = intval(AF_2K3N_mat);
end    
AF_2K3N_mat(Ind_notKin2K,:) = 1i * AntiDer_2K_tail * [F_2K2N_mat(Ind_notKin2K,:), zeros(l_notKin2K,N-1)] * transpose(MUpsilon_3N); % tail part

%ADF(X)
DF_32K2N_tens = DF_Lorenz_tens(X_3K2N, para_2N); %% /!\ We only need the columns up to |k| = 2K-2, but these columns have non zero entries up to the row 3K-3
DF_32K2N_tens = DF_32K2N_tens(:,Ind_2Kin3K,:);
DF_12K2N_mat = reshape( permute( DF_32K2N_tens(Ind_Kin3K,:,:), [3,1,2] ), [N2*(6*K-2),6*K2-2]); % finite part of DF which will be multiplied by A later on
ADF_32K2N_tens = zeros(6*K3-2,6*K2-2,N2);
if exist('intval','file') && isintval(X(1))
    ADF_32K2N_tens = intval(ADF_32K2N_tens);
end    
ADF_32K2N_tens(Ind_notKin3K,:,:) = permute( reshape( MUpsilon_2N * reshape( permute( DF_32K2N_tens(Ind_notKin3K,:,:), [3,1,2] ), [N2,(6*K2-2)*l_notKin3K] ), [N2,l_notKin3K,6*K2-2]), [2,3,1]);
ADF_32K2N_tens(Ind_notKin3K,:,:) = 1i * reshape( AntiDer_3K_tail * reshape( ADF_32K2N_tens(Ind_notKin3K,:,:), [l_notKin3K,N2*(6*K2-2)] ), [l_notKin3K,6*K2-2,N2] ); % tail part

for ind = 1:6*K-2 % multiplication by a for the finite part of F and DF, row by row (in k) for a   
    %AF(X)
    A_k3N = convo_PC_mat( squeeze(a_K3N(ind,:,:)), alpha_mat_3N, 'column');
    AF_2K3N_mat(Ind_Kin2K(ind),:) = A_k3N * F_K3N_vect;
    %ADF(X)
    A_k2N = convo_PC_mat( squeeze(a_K3N(ind,:,1:N2)), alpha_mat_2N, 'column');
    ADF_32K2N_tens(Ind_Kin3K(ind),:,:) = permute( reshape( A_k2N * DF_12K2N_mat, [N2,1,6*K2-2] ), [2,3,1] );  
end

Y = weights_2K3N * abs(tens2vect(AF_2K3N_mat));
fprintf('\nY = %g\n',i2f(Y))


B_tens = I_tens(:,Ind_2Kin3K,:) - ADF_32K2N_tens;
Z1_finite = norm_op_tens(B_tens, weights_2N', weights_3K, weights_2K);
fprintf('\nZ1_finite = %g\n',i2f(Z1_finite))


%% Z1_tail
x_K2N = [x zeros(2*K-1,N-1)];
y_K2N = [y zeros(2*K-1,N-1)];
z_K2N = [z zeros(2*K-1,N-1)];
MUpsilon_2N_diag = convo_PC_diag_mat([Upsilon; zeros(N-1,1)], 2*K-1, alpha_mat_2N);
e_2N = [e; zeros(N-1,1)];
Omega_2N = [transpose(Omega); zeros(N-1,1)];

norm_emUpsilonOmega = weights_2N * abs(e_2N - MUpsilon_2N*Omega_2N);
norm_Upsilon = weights_N * abs(Upsilon);
norm_xUpsilon = weights_singleK2N * abs( MUpsilon_2N_diag * tens2vect(x_K2N) );
norm_yUpsilon = weights_singleK2N * abs( MUpsilon_2N_diag * tens2vect(y_K2N) );
norm_zmrhoUpsilon = weights_singleK2N * abs( MUpsilon_2N_diag * tens2vect( z_K2N - [zeros(K-1,N2);rho_2N;zeros(K-1,N2)])  );
Z1_tail = norm_emUpsilonOmega + ...
          1/K * max([ sigma*norm_Upsilon + norm_zmrhoUpsilon + norm_yUpsilon;
                      (sigma + 1)*norm_Upsilon + norm_xUpsilon;
                      norm_xUpsilon + beta*norm_Upsilon]);
fprintf('\nZ1_tail = %g\n',i2f(Z1_tail))
            
Z1 = max(Z1_finite,Z1_tail);
            
%% Z2
Der_K = spdiags( [0 repmat(abs(-K+1:K-1), [1,3])]', 0, 6*K-2, 6*K-2);
ak = permute(reshape( reshape(permute(a,[1,3,2]),[N*(6*K-2),6*K-2]) * Der_K, [6*K-2,N,6*K-2]),[1,3,2]);
Z2 = max([ norm_op_tens(a, weights_N', weights_K);
           norm_op_tens(ak, weights_N', weights_K);
           norm_Upsilon; ]);
       
fprintf('\nZ2 = %g\n',i2f(Z2))

%% Cheking that we have a contraction
[rmin, rmax] = iscontraction(Y, Z1, Z2);

