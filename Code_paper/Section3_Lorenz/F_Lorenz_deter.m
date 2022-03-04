function FF = F_Lorenz_deter(X,para)

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

xy = convo_exp(x,y);
xz = convo_exp(x,z);
MK = diag(-K+1:K-1);

G = 1i* ( ((-K+1:K-1).*xphase')*x + ((-K+1:K-1).*yphase')*y + ((-K+1:K-1).*zphase')*z );
Fx = -1i*MK*Omega*x + sigma*(y-x);
Fy = -1i*MK*Omega*y + rho*x - y - xz;
Fz = -1i*MK*Omega*z -beta*z + xy;

FF = [G; Fx; Fy; Fz];


