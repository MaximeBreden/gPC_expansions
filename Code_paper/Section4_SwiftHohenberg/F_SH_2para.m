function F = F_SH_2para(u, para)

L = para.L;
alpha_mat = para.alpha_mat;

[K,N1,N2] = size(i2f(u));

Mu = convo_cos_PC_mat(u,alpha_mat);
uuu_vect = Mu*(Mu*tens2vect(u));

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

F = vect2tens((Mrho - Lambda)*tens2vect(u) - Mbeta*uuu_vect, [K,N1,N2]);


