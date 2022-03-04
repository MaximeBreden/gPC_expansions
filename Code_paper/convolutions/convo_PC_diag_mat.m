function Mb_diag = convo_PC_diag_mat(b,nb,alpha_mat)

N = my_numel(b);
Mb = convo_PC_mat(b,alpha_mat);
Mb_mult = repmat(Mb,1,nb);                                   
Mb_cell = mat2cell(Mb_mult,N,repmat(N,1,nb));  
Mb_diag = blkdiag(Mb_cell{:});  