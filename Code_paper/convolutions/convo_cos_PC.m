function c = convo_cos_PC(a,b,alpha_mat)

Ma = convo_cos_PC_mat(a, alpha_mat);
c = vect2tens( Ma*tens2vect(b), size(b) );

% Ma = convo_cos_PC_mat(a, alpha_mat);
% bb = reshape(transpose(b), numel(b), 1);
% c = transpose( reshape( Ma*bb, flip(size(b)) ) );