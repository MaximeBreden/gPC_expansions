function Mu = convo_cos_PC_mat(u,alpha_mat)

[K,N1,N2] = size(i2f(u,'mid'));
H_K = hankel((1:K)');
T_K = toeplitz((1:K)');
T_K(:,1) = 0;

u = cat(1, zeros(1,N1,N2), u);
extra_dim = 1;
Mu = convo_PC_mat(u,alpha_mat,extra_dim,'column');
N = N1*N2;

if exist('intval','file') && isintval(u(1))
    % The mat2cell/cell2mat conversion can be made to work directly with
    % intvals, but it is painfully slow, so we proceed as follows instead
    Mu_mid = mat2cell(mid(Mu), N, N*ones(1,K+1));
    MuT_mid = cell2mat(Mu_mid(T_K+1));
    MuH_mid = cell2mat(Mu_mid(H_K+1));
    Mu_rad = mat2cell(rad(Mu), N, N*ones(1,K+1));
    MuT_rad = cell2mat(Mu_rad(T_K+1));
    MuH_rad = cell2mat(Mu_rad(H_K+1));
    MuT = midrad(MuT_mid,MuT_rad);
    MuH = midrad(MuH_mid,MuH_rad);
    Mu = MuT + MuH;
else
    Mu = mat2cell(Mu, N, N*ones(1,K+1));
    Mu = cell2mat(Mu(H_K+1)) + cell2mat(Mu(T_K+1));
end   


% The above code generalized what is below to potentially 2D gPC extensions
%
% [K,N] = size(u);
% H_K = hankel((1:K)');
% T_K = toeplitz((1:K)');
% T_K(:,1) = 0;
% 
% Mu = convo_PC_mat([zeros(1,N);u],alpha_mat,'column');
% Mu = mat2cell(Mu, N, N*ones(1,K+1));
% Mu = my_cell2mat(Mu(H_K+1)) + my_cell2mat(Mu(T_K+1));

