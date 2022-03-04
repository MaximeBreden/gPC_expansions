function DF = DF_SH_2para(u, para)

L = para.L;
alpha_mat = para.alpha_mat;

[K,N1,N2] = size(i2f(u));

Muu = convo_cos_PC_mat(convo_cos_PC(u,u,alpha_mat), alpha_mat);
Lambda = (1-((0:K-1)'*pi/L).^2).^2;
Lambda = spdiags(repelem(Lambda,N1*N2,1), 0, K*N1*N2, K*N1*N2);

rho = zeros(K,N1,N2);
if exist('intval','file') && isintval(u(1))
    rho = intval(rho);
end 
rho(1,:,1) = para.rho;
Mrho = convo_cos_PC_mat(rho,alpha_mat);
beta = zeros(K,N1,N2);
if exist('intval','file') && isintval(u(1))
    beta = intval(beta);
end 
beta(1,1,:) = para.beta;
Mbeta = convo_cos_PC_mat(beta,alpha_mat);

DF = Mrho - Lambda - 3*Mbeta*Muu;