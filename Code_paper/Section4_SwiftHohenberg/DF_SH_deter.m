function DF = DF_SH_deter(u, para)

L = para.L;
rho = para.rho;

K = length(u);

Muu = convo_cos_mat(convo_cos(u,u));
Lap = diag( -((0:K-1)*pi/L).^2 );
I = eye(K);

DF = (rho*I-(I+Lap)^2) - 3*Muu;