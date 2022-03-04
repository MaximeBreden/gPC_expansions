function Mu = convo_exp_mat(u)

K = (length(u)+1)/2;
Mu = toeplitz( transpose([u(K:2*K-1);zeros(K-1,1)]), [u(K:-1:1);zeros(K-1,1)] );
