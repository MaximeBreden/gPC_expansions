function Mb = convo_PC_mat_1D(b,alpha_mat,output_format)

% If b is a vector, the output is simply the associated generalized 
% convolution matrix (for the linearization coefficients in alpha_mat).

% If b is a matrix with rows (or columns) of PC coefficients, return the
% collection of convolution matrices associated to each row (or column).
% Takes coefficients row by row by default, but switches to column by
% column if the size of the input is only compatible this way. The output
% convolution matrices are concatenated "row-wise" or "column-wise"
% depending on the value of the input parameter format ("column-wise" is 
% selected by default).


N = size(alpha_mat,2);
[L,C] = size(b);
if C==N && L~=N
    b = transpose(b);
    C = L;
end

if nargin==2 
    output_format = 'column';
end

Mb = alpha_mat*b;

switch output_format
    case 'column'
        Mb = reshape(Mb,[N,N*C]);
    case 'row'
        Mb = reshape( permute( reshape(Mb,[N,N,C]), [1,3,2]), [N*C,N] );
    otherwise
        error('This format is not valid')
end




