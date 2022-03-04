function DF = DF_SH_2para_tens(u, para)

L = para.L;
alpha_mat = para.alpha_mat;

[K,N1,N2] = size(i2f(u));

Lambda = diag((1-((0:K-1)'*pi/L).^2).^2);
Mlambda = zeros(K,K,N1,N2);
Mlambda(:,:,1,1) = Lambda;

rho = zeros(K,N1,N2);
if exist('intval','file') && isintval(u(1))
    rho = intval(rho);
end 
rho(1,:,1) = para.rho;
Mrho = convo_cos_PC_tens(rho);
beta = zeros(K,N1,N2);
if exist('intval','file') && isintval(u(1))
    beta = intval(beta);
end 
beta(1,1,:) = para.beta;
Mbetauu = convo_cos_PC_tens(convo_cos_PC(beta,convo_cos_PC(u,u,alpha_mat),alpha_mat));

DF = Mrho - Mlambda - 3*Mbetauu;