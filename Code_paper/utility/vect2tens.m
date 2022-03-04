function v = vect2tens(u,sz)

% Inverse function of tens2vect.

v = permute( reshape(u,flip(sz)), length(sz):-1:1 );
