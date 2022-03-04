function DF = DF_SH(u, para)

L = para.L;
beta = para.beta;
rho = para.rho;
alpha_mat = para.alpha_mat;

[K,N] = size(u);

Muu = convo_cos_PC_mat(convo_cos_PC(u,u,alpha_mat), alpha_mat);
Mrho_diag = convo_PC_diag_mat(transpose(rho), K, alpha_mat);
Lambda = (1-((0:K-1)*pi/L).^2).^2;
Lambda_N = spdiags( repelem(Lambda, N)', 0, K*N, K*N );

DF = Mrho_diag - Lambda_N - 3*beta*Muu;