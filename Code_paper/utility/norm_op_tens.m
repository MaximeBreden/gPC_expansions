function val = norm_op_tens(a, weights_N, weights_K_row, weights_K_column)

% Computes a matrix Mat with entries equal to the norm of each element 
% a(i,j,{:}), using the weights in weights_N, and then computes the 
% operator norm of Mat, using the weights in weights_K_row and
% weights_K_column

if nargin < 4
    weights_K_column = weights_K_row;
end

l_r = length(weights_K_row);
l_c = length(weights_K_column);

weights_N_tens = superkron(ones(l_r,1,1,1),ones(1,l_c,1,1),reshape(weights_N,[1,1,size(weights_N)]));
sum_dim = 3:my_ndims(a);
if isempty(sum_dim) %degenerate case where the extension in N is trivial (N=1)
    Mat = abs(a) .* weights_N_tens;
else
    Mat = sum(abs(a) .* weights_N_tens, sum_dim);
end
val = max( (weights_K_row * Mat) ./ weights_K_column );


%%% Perhaps slightly easier to read code, for the particular case where each
%%% a(i,j,{:}) is only one-dimensional
% if nargin < 4
%     weights_K_column = weights_K_row;
% end
% if size(weights_N,1) == 1
%     weights_N = transpose(weights_N);
% end
% 
% l_r = length(weights_K_row);
% l_c = length(weights_K_column);
% 
% N = length(weights_N);
% weights_N_tens = permute( reshape( repmat( weights_N, [l_r,l_c] ), [N,l_r,l_c] ), [2,3,1] );
% Mat = sum( abs(a) .* weights_N_tens, 3);
% val = max( (weights_K_row * Mat) ./ weights_K_column );