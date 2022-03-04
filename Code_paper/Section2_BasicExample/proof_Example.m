function [rmin, rmax, a] = proof_Example(x, para, a)

%The last argument is optional, and will be computed inside this function
%if not given as input.

if exist('intval','file') && isintval(x(1))
    fprintf('\nRigorous validation with interval arithmetic ...\n')
else
    fprintf('\nPrevalidation without interval arithmetic ...\n')
end

%% Initialization
g = para.g;
eta = para.eta;
PC_type = para.PC_type;
mu = para.mu;
N = length(x);
e = zeros(N, 1);
e(1) = 1;

%% Building A
if nargin < 3
    a = DF_Example(x, para) \ e; 
end

%% Zero padding to remove truncation errors
N2 = 2*N-1;
N3 = 3*N-2;

x_2N = [x; zeros(N-1,1)];
x_3N = [x; zeros(2*(N-1),1)];
e_2N = [e; zeros(N-1,1)];
a_2N = [a; zeros(N-1,1)];
a_3N = [a; zeros(2*(N-1),1)];
g_3N = [g; zeros(2*(N-1),1)];

weights_N = eta.^(0:N-1);
weights_2N = eta.^(0:N2-1);
weights_3N = eta.^(0:N3-1);

alpha_tens_3N = Tens_for_prod(N3,PC_type,mu);
alpha_mat_3N = reshape(alpha_tens_3N,[N3^2,N3]);
alpha_mat_2N = reshape(alpha_tens_3N(1:N2,1:N2,1:N2),[N2^2,N2]);

para_3N = para;
para_3N.g = g_3N;
para_3N.alpha_mat = alpha_mat_3N;

%% Y
F_3N = F_Example(x_3N, para_3N);
AF_3N = convo_PC(a_3N, F_3N, alpha_mat_3N);
Y = weights_3N * abs(AF_3N)

%% Z1
b_2N = e_2N - convo_PC(a_2N, 2*x_2N, alpha_mat_2N);
Z1 = weights_2N * abs(b_2N)
       
%% Z2
Z2 = weights_N * abs(2*a)

%% Cheking that we have a contraction
[rmin, rmax] = iscontraction(Y, Z1, Z2);
