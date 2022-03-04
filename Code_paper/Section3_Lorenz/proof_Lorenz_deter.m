function [rmin, rmax] = proof_Lorenz_deter(X,para)

sigma = para.sigma;
rho = para.rho;
beta = para.beta;
nu = para.nu;

K = (length(X)+2)/6;
Omega = X(1);
x = X(2:2*K);
y = X(2*K+1:4*K-1);
z = X(4*K:6*K-2);

K2 = 2*K-1;
K3 = 3*K-2;
X_2K = resize_deter(X,K2);
X_3K = resize_deter(X,K3);
weights_K = [1 nu.^repmat(abs(-K+1:K-1), [1,3])];
weights_2K = [1 nu.^repmat(abs(-K2+1:K2-1), [1,3])];
weights_3K = [1 nu.^repmat(abs(-K3+1:K3-1), [1,3])];

%%% Building A
DF = DF_Lorenz_deter(X,para);
A_K = inv(DF);

Ind_K_single = -K+1:K-1;
Ind_K = [1, K3+1+Ind_K_single, 3*K3+Ind_K_single, 5*K3-1+Ind_K_single];
A_3K = diag(1./[1 repmat(-K3+1:K3-1, [1,3])]) / (-1i*Omega);
A_3K(Ind_K,Ind_K) = A_K;

%%% Y
F_2K = F_Lorenz_deter(X_2K,para);
Ind_2K_single = -K2+1:K2-1;
Ind_2K = [1, K3+1+Ind_2K_single, 3*K3+Ind_2K_single, 5*K3-1+Ind_2K_single];
AF_2K = A_3K(Ind_2K,Ind_2K) * F_2K;
Y = weights_2K * abs(AF_2K)


%%% Z1
I_3K = eye(6*K3-2);
DF_3K = DF_Lorenz_deter(X_3K,para);
B = I_3K(:,Ind_2K)- A_3K*DF_3K(:,Ind_2K);
Z1_finite = max( (weights_3K * abs(B)) ./ weights_2K )

weights_single = nu.^abs(-K+1:K-1);
e = zeros(2*K-1,1);
e(K) = 1;
Z1_tail = max([ sigma + weights_single*abs(z-rho*e) + weights_single*abs(y);
                sigma + 1 + weights_single*abs(x);
                weights_single*abs(x) + beta]) / (Omega*K)
            
Z1 = max(Z1_finite,Z1_tail);
            
%%% Z2
Der = [0 repmat(abs(-K+1:K-1), [1,3])];
Z2 = max([ max( (weights_K * abs(A_K*diag(Der))) ./ weights_K );
           1/Omega;      
           max( (weights_K * abs(A_K)) ./ weights_K ) ])

%% Cheking that we have a contraction
[rmin, rmax] = iscontraction(Y, Z1, Z2);

