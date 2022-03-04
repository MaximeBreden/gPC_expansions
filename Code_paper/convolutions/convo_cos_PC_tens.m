function Mu = convo_cos_PC_tens(u)

[K,N1,N2] = size(i2f(u,'mid'));
H_K = hankel((1:K)');
T_K = toeplitz((1:K)');
T_K(:,1) = 0;

if exist('intval','file') && isintval(u(1))
    % The mat2cell/cell2mat conversion can be made to work directly with
    % intvals, but it is painfully slow, so we proceed as follows instead
    Mu_mid = mat2cell( reshape( cat(1,zeros(1,N1,N2),mid(u)), [1,K+1,N1,N2] ), 1, ones(1,K+1), N1, N2 );
    MuT_mid = cell2mat(Mu_mid(T_K+1));
    MuH_mid = cell2mat(Mu_mid(H_K+1));
    Mu_rad = mat2cell( reshape( cat(1,zeros(1,N1,N2),rad(u)), [1,K+1,N1,N2] ), 1, ones(1,K+1), N1, N2 );
    MuT_rad = cell2mat(Mu_rad(T_K+1));
    MuH_rad = cell2mat(Mu_rad(H_K+1));
    MuT = midrad(MuT_mid,MuT_rad);
    MuH = midrad(MuH_mid,MuH_rad);
    Mu = MuT + MuH;
else
    Mu = mat2cell( reshape( cat(1,zeros(1,N1,N2),u), [1,K+1,N1,N2] ), 1, ones(1,K+1), N1, N2 );
    Mu = cell2mat(Mu(H_K+1)) + cell2mat(Mu(T_K+1));
end   



% [K,N] = size(u);
% H_K = hankel((1:K)');
% T_K = toeplitz((1:K)');
% T_K(:,1) = 0;
% 
% Mu = mat2cell( reshape( [zeros(1,N);u], [1,K+1,N] ), 1, ones(1,K+1), N );
% Mu = cell2mat(Mu(H_K+1)) + cell2mat(Mu(T_K+1));