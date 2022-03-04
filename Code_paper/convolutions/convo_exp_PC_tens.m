function Mu = convo_exp_PC_tens(u)

[K,N] = size(u);
K = (K+1)/2;
T_K = toeplitz([(K:2*K-1) zeros(1,K-1)]',[(K:-1:1) zeros(1,K-1)]);

if exist('intval','file') && isintval(u(1))
    % The mat2cell/cell2mat conversion can be made to work directly with
    % intvals, but it is painfully slow, so we proceed as follows instead
    Mu_mid = mat2cell( reshape( [zeros(1,N);mid(u)], [1,2*K,N] ), 1, ones(1,2*K), N );
    Mu_mid = cell2mat(Mu_mid(T_K+1));
    Mu_rad = mat2cell( reshape( [zeros(1,N);rad(u)], [1,2*K,N] ), 1, ones(1,2*K), N );
    Mu_rad = cell2mat(Mu_rad(T_K+1));
    Mu = midrad(Mu_mid,Mu_rad);
else
    Mu = mat2cell( reshape( [zeros(1,N);u], [1,2*K,N] ), 1, ones(1,2*K), N );
    Mu = cell2mat(Mu(T_K+1));
end   

