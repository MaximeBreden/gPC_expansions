function DFF = DF_Lorenz_deter(X,para)

sigma = para.sigma;
rho = para.rho;
beta = para.beta;
Xphase = para.Xphase;

K = (length(X)+2)/6;
Omega = X(1);
x = X(2:2*K);
y = X(2*K+1:4*K-1);
z = X(4*K:6*K-2);

Kphase = (length(Xphase)+3)/6;
xphase = [zeros(K-Kphase,1); Xphase(1:2*Kphase-1); zeros(K-Kphase,1)];
yphase = [zeros(K-Kphase,1); Xphase(2*Kphase:4*Kphase-2); zeros(K-Kphase,1)];
zphase = [zeros(K-Kphase,1); Xphase(4*Kphase-1:6*Kphase-3); zeros(K-Kphase,1)];

I = eye(2*K-1);
MK = diag(-K+1:K-1);
Mx = convo_exp_mat(x);
My = convo_exp_mat(y);
Mz = convo_exp_mat(z);

DGDOmega = 0;
DFxDOmega = -1i*MK*x;
DFyDOmega = -1i*MK*y;
DFzDOmega = -1i*MK*z;

DGDx = 1i*(-K+1:K-1).*xphase';
DGDy = 1i*(-K+1:K-1).*yphase';
DGDz = 1i*(-K+1:K-1).*zphase';

DFxDx = -1i*Omega*MK - sigma*I;
DFxDy = sigma*I;
DFxDz = zeros(2*K-1);

DFyDx = rho*I - Mz;
DFyDy = -1i*Omega*MK - I;
DFyDz = -Mx;

DFzDx = My;
DFzDy = Mx;
DFzDz = -1i*Omega*MK - beta*I;

DFF = [DGDOmega DGDx DGDy DGDz
       DFxDOmega DFxDx DFxDy DFxDz
       DFyDOmega DFyDx DFyDy DFyDz
       DFzDOmega DFzDx DFzDy DFzDz];