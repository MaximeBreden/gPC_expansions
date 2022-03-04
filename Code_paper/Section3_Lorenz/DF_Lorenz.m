function DF = DF_Lorenz(X, para)

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

I = eye((2*K-1)*N);
MK = spdiags( repelem(-K+1:K-1, N)', 0, (2*K-1)*N, (2*K-1)*N );
Mx = convo_exp_PC_mat(x,alpha_mat);
My = convo_exp_PC_mat(y,alpha_mat);
Mz = convo_exp_PC_mat(z,alpha_mat);
 
MOmega_diag = convo_PC_diag_mat(Omega,2*K-1,alpha_mat);
Mrho_diag = convo_PC_diag_mat(rho,2*K-1,alpha_mat);

DG = 1i * convo_PC_mat( spdiags([0; repmat((-K+1:K-1)',3,1)], 0, 6*K-2, 6*K-2) * conj(Xphase) , alpha_mat, 1, 'column');
DFDOmega = -1i * convo_PC_mat( spdiags(repmat((-K+1:K-1)',3,1), 0, 6*K-3, 6*K-3) * X(2:end,:) , alpha_mat, 1, 'row');

DFxDx = -1i*MK*MOmega_diag - sigma*I;
DFxDy = sigma*I;
DFxDz = zeros((2*K-1)*N);

DFyDx = Mrho_diag - Mz;
DFyDy = -1i*MK*MOmega_diag - I;
DFyDz = -Mx;

DFzDx = My;
DFzDy = Mx;
DFzDz = -1i*MK*MOmega_diag - beta*I;


DF = [DG
      DFDOmega [DFxDx DFxDy DFxDz
                DFyDx DFyDy DFyDz
                DFzDx DFzDy DFzDz]];