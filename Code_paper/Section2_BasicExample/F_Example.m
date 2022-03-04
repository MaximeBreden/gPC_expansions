function F = F_Example(x, para)

g = para.g;
alpha_mat = para.alpha_mat;

F = convo_PC(x,x,alpha_mat) - g;