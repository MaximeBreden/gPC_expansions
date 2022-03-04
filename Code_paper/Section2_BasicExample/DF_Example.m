function DF = DF_Example(x, para)

alpha_mat = para.alpha_mat;

DF = 2*convo_PC_mat(x,alpha_mat);