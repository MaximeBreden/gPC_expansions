function F = F_Lorenz(X, para)

sigma = para.sigma;
rho = para.rho;
beta = para.beta;
Xphase = para.Xphase;
alpha_mat = para.alpha_mat;

[K,N] = size(X);
K = (K+2)/6;

Omega = X(1,:);
x = X(2:2*K,:);
y = X(2*K+1:4*K-1,:);
z = X(4*K:6*K-2,:);

Xphase = [zeros(1,size(Xphase,2)); Xphase];
Xphase = resize(Xphase, K, N);

xy = convo_exp_PC(x,y,alpha_mat);
xz = convo_exp_PC(x,z,alpha_mat);
xrho = x * transpose( convo_PC_mat(rho,alpha_mat) );
MOmega = transpose( convo_PC_mat(Omega,alpha_mat) );
MK = diag(-K+1:K-1);


G = 1i * transpose( convo_PC_mat( spdiags([0; repmat((-K+1:K-1)',3,1)], 0, 6*K-2, 6*K-2) * conj(Xphase) , alpha_mat, 1, 'column') * tens2vect(X) );
Fx = -1i*MK*x*MOmega + sigma*(y-x);
Fy = -1i*MK*y*MOmega + xrho - y - xz;
Fz = -1i*MK*z*MOmega -beta*z + xy;

F = [G; Fx; Fy; Fz];

if isfield(para,'reshapeX')
    F = para.reshapeX(F);
end

