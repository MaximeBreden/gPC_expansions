function DF = DF_SH_tens(u, para)

L = para.L;
beta = para.beta;
rho = para.rho;
alpha_mat = para.alpha_mat;

[K,N] = size(u);

Muu = convo_cos_PC_tens( convo_cos_PC(u, u, alpha_mat) );
Mrho = convo_cos_PC_tens( [rho; zeros(K-1,N)] );
Mlambda = zeros(K, K, 1);
Mlambda(:) = diag((1-((0:K-1)*pi/L).^2).^2);
Mlambda = cat(3, Mlambda, zeros(K, K, N-1));

DF = Mrho - Mlambda - 3*beta*Muu;