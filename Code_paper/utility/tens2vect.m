function v = tens2vect(u)

% Converts a tensor of Fourier x gPC coefficients into a vector, with the 
% suitable ordering for matrix/vector multiplications. The gPC expansions
% itself can have several dimensions.

v = reshape( permute( u, my_ndims(u):-1:1 ), my_numel(u), 1 );
