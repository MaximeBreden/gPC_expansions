function K = my_kron(A,B)

% Kronecker product, compatible with matrices of intvals

[ma,na] = size(A);
[mb,nb] = size(B);
[ia,ib] = meshgrid(1:ma,1:mb);
[ja,jb] = meshgrid(1:na,1:nb);
K = A(ia,ja).*B(ib,jb);