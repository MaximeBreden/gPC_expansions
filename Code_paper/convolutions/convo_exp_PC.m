function c = convo_exp_PC(a,b,alpha_mat)

Ma = convo_exp_PC_mat(a, alpha_mat);
c = vect2tens( Ma*tens2vect(b), size(b) );