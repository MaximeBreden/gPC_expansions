function F = F_SH(u, para)

L = para.L;
beta = para.beta;
rho = para.rho;
alpha_mat = para.alpha_mat;

[K,~] = size(u);

Mu = convo_cos_PC_mat(u,alpha_mat);
uuu = vect2tens( Mu*(Mu*tens2vect(u)), size(u) );
% uuu = convo_cos_PC(u, convo_cos_PC(u,u,alpha_mat), alpha_mat);

Lambda = diag((1-((0:K-1)*pi/L).^2).^2);
Mrho = transpose( convo_PC_mat(rho,alpha_mat,1) );

F = u*Mrho - Lambda*u - beta*uuu;

