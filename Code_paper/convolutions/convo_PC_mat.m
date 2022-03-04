function Mb = convo_PC_mat(b, alpha_mat, extra_dim, output_format)

% Constructs the generalized convolution matrix for b associated to the 
% linearization coefficients in alpha_mat. 
%
% b can be a vector, or a higher dimensional tensor, and alpha_mat can be
% a matrix or a cell array of matrices. If extra_dim > 0, then we do not
% consider generalized convolutions along the extra_dim first dimensions in
% b. For instance if b is a tensor containing Fourier x gPC coefficients,
% by putting extra_dim = 1 we produce generalized convolution matrices for
% each Fourier coefficient in b. In such a case, output_format selects in
% which way these matrices are concatenated.
%
% The generic function for this is convo_PC_mat_anyD, but we use better 
% optimized code depending on the cases.

if nargin < 4
    output_format = 'column';
    if nargin < 3
        extra_dim = 0;
    end
end

if not(iscell(alpha_mat))
    Mb = convo_PC_mat_1D(b,alpha_mat,output_format);
elseif length(alpha_mat) == 1
    Mb = convo_PC_mat_1D(b,alpha_mat{:},output_format);
elseif length(alpha_mat) == 2 
    if alpha_mat{2} == 1
        Mb = convo_PC_mat_1D(b,alpha_mat{1},output_format);
    elseif ismatrix(b) && extra_dim == 0
        Mb = convo_PC_mat_2D(b,alpha_mat{1},alpha_mat{2});
    else
        Mb = convo_PC_mat_anyD_opt(b,alpha_mat,extra_dim,output_format);
    end
else
    Mb = convo_PC_mat_anyD_opt(b,alpha_mat,extra_dim,output_format);
end

if nnz(i2f(Mb)) / numel(i2f(Mb)) < 0.1
    Mb = sparse(Mb);
end