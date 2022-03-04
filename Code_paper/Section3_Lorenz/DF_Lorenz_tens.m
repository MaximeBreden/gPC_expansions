function DF = DF_Lorenz_tens(X, para)

sigma = para.sigma;
rho = para.rho;
beta = para.beta;
Xphase = para.Xphase;

[K,N] = size(X);
K = (K+2)/6;

Omega = X(1,:);
x = X(2:2*K,:);
y = X(2*K+1:4*K-1,:);
z = X(4*K:6*K-2,:);

Xphase = [zeros(1,size(Xphase,2)); Xphase];
Xphase = resize(Xphase, K, N);

e = [1 zeros(1,N-1)];

I = zeros(2*K-1, 2*K-1, 1);
I(:) = eye(2*K-1);
I = cat(3, I, zeros(2*K-1, 2*K-1, N-1));

Mx = convo_exp_PC_tens(x);
My = convo_exp_PC_tens(y);
rhomz = [zeros(K-1,N);rho;zeros(K-1,N)]-z;
Mrhomz = convo_exp_PC_tens(rhomz);


%%
DG = 1i * reshape( spdiags([0; repmat((-K+1:K-1)',3,1)], 0, 6*K-2, 6*K-2) * conj(Xphase), [1,6*K-2,N]);
DFDOmega = -1i * reshape( spdiags(repmat((-K+1:K-1)',3,1), 0, 6*K-3, 6*K-3) * X(2:end,:), [6*K-3,1,N]);

DFxDx = zeros(2*K-1, 2*K-1, N);
DFyDy = zeros(2*K-1, 2*K-1, N);
DFzDz = zeros(2*K-1, 2*K-1, N);
if exist('intval','file') && isintval(X(1))
    DFxDx = intval(DFxDx);
    DFyDy = intval(DFyDy);
    DFzDz = intval(DFzDz);
end
for k = 1:2*K-1
    DFxDx(k,k,:) = -1i*(k-K)*Omega - sigma*e;
    DFyDy(k,k,:) = -1i*(k-K)*Omega - e;
    DFzDz(k,k,:) = -1i*(k-K)*Omega - beta*e;
end
DFxDy = sigma*I;
DFxDz = zeros(2*K-1, 2*K-1, N);

DFyDx = Mrhomz;
DFyDz = -Mx;

DFzDx = My;
DFzDy = Mx;

DF = cat(1, DG,...
            cat(2, DFDOmega, cat(1, cat(2, DFxDx, DFxDy, DFxDz),...
                                    cat(2, DFyDx, DFyDy, DFyDz),...
                                    cat(2, DFzDx, DFzDy, DFzDz) ) ) );
                                
% % DF = [DG
% %       DFDOmega [DFxDx DFxDy DFxDz
% %                 DFyDx DFyDy DFyDz
% %                 DFzDx DFzDy DFzDz]];