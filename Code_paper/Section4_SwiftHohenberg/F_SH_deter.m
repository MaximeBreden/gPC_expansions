function F = F_SH_deter(u, para)

L = para.L;
rho = para.rho;

K = length(u);

Mu = convo_cos_mat(u);
uuu = Mu*Mu*u;
Lap = diag( -((0:K-1)*pi/L).^2 );
I = eye(K);

F = (rho*I-(I+Lap)^2)*u -uuu;


