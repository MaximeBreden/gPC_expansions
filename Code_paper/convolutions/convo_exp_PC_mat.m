function Mu = convo_exp_PC_mat(u,alpha_mat)

[K,N] = size(u);
K = (K+1)/2;
T_K = toeplitz([(K:2*K-1) zeros(1,K-1)]',[(K:-1:1) zeros(1,K-1)]);

Mu = convo_PC_mat([zeros(1,N);u],alpha_mat,'column');
if exist('intval','file') && isintval(u(1))
    % The mat2cell/cell2mat conversion can be made to work directly with
    % intvals, but it is painfully slow, so we proceed as follows instead
    Mu_mid = mat2cell(mid(Mu), N, N*ones(1,2*K));
    Mu_mid = cell2mat(Mu_mid(T_K+1));
    Mu_rad = mat2cell(rad(Mu), N, N*ones(1,2*K));
    Mu_rad = cell2mat(Mu_rad(T_K+1));
    Mu = midrad(Mu_mid,Mu_rad);
else
    Mu = mat2cell(Mu, N, N*ones(1,2*K));
    Mu = cell2mat(Mu(T_K+1));
end   
