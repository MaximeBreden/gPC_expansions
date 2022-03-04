function [rmin, rmax] = proof_SH_deter(u,para,rstar)

L = para.L;
rho = para.rho;
nu = para.nu;

K = length(u);
weights = [1 2*ones(1,K-1)] .* nu.^(0:K-1);
K3 = 3*K-2;
K5 = 5*K-4;
weights_K5 = [1 2*ones(1,K5-1)] .* nu.^(0:K5-1);
weights_K3 = weights_K5(1:K3);
u_K3 = [u;zeros(2*(K-1),1)];
u_K5 = [u;zeros(4*(K-1),1)];

%%% Building A (lazy)
DF = DF_SH_deter(u,para);
A = inv(DF);

Lambda_K5 = - (1-((0:K5-1)*pi/L).^2).^2;
A_K5 = diag(1./Lambda_K5);
A_K5(1:K,1:K) = A;

%%% Y
F_K3 = F_SH_deter(u_K3,para);
AF_K3 = A_K5(1:K3,1:K3) * F_K3;
Y = weights_K3 * abs(AF_K3)


%%% Z1
DF_K5 = DF_SH_deter(u_K5,para);
I = eye(K5);
B = I(:,1:K3) - A_K5*DF_K5(:,1:K3);
Z1_finite = max( (weights_K5 * abs(B)) ./ weights_K3 )

if (K*pi/L)^2 >= 1
    min_lambda_tail = (1-(K*pi/L)^2)^2;
else
    warning('K may not be large enough')
    K_lim = ceil( L/pi );
    min_lambda_tail = min( (1-((K:K_lim)*pi/L).^2).^2 );
end
uu_K3 = convo_cos(u_K3,u_K3);
e_K3 = zeros(K3,1);
e_K3(1) = 1;
Z1_tail = weights_K3 * abs(rho*e_K3-3*uu_K3) / min_lambda_tail
            
Z1 = max(Z1_finite,Z1_tail);
            
%%% Z2
norm_A_finite = max((weights*abs(A)./weights));
Z2 = max(norm_A_finite, 1/min_lambda_tail ) * 6 * ( weights*abs(u) + rstar )

[rmin, rmax] = iscontraction(Y, Z1, Z2, rstar);
