function Mb = convo_PC_mat_2D(b,alpha1_mat,alpha2_mat)

% Optimized 2D generalized convolution matrix, works only if the dimension
% of b is equal to 2.

[N1,N2] = size(b);
Mb2 = convo_PC_mat_1D(transpose(b),alpha2_mat,'row');
Mb = my_kron(alpha1_mat,speye(N2)) * Mb2;

Mb = mat2cell(Mb, N2*ones(1,N1^2), N2); % a column of N1^2 N2xN2 blocks
Mb = reshape(Mb, N1, N1);
Mb = cell2mat(Mb);

