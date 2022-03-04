function Tens = Tens_for_prod(N,PC_type,arg)

%Construct the 3D tensor Tens such that, whenever b is a column vector of
%size N, Mb=reshape(reshape(Tens,[N^2,N])*b,[N,N]) is the PC product matrix
%of b (s.t. for any vector a of size N, prod_PC(a,b)=Mb*a).
%The last argument is optional. If it contains a number of intval type, the
%linearization coefficients are computed using interval arithmetic. For
%Gegenbauer polynomial, this last argument is mandatory and corresponds to
%the parameter mu.

if nargin == 2
    arg = [];
end

%% We first check whether the required linearization coefficients were already precomputed
if exist('intval','file') && isintval(arg)
    str = ['linearization_coefficients/LinearizationCoeffs_',PC_type,'_Intval.mat'];
else
    str = ['linearization_coefficients/LinearizationCoeffs_',PC_type,'.mat'];
end
if exist(str,'file')==2
    load(str,'Tens')
    if size(Tens,1) >= N
        Tens = Tens(1:N,1:N,1:N);
        return
    end
end


%% If needed, we compute and store the required linearization coefficients
fprintf("\nThis computation requires linearization coefficients that have not been precomputed yet, we do so now ...\n")
[A,B,C] = meshgrid(0:N-1,0:N-1,0:N-1);
C = min(A,min(B,C));
switch PC_type
    case 'Legendre'
        Tens_alpha = alpha_L(A,B,C,arg);
    case 'Chebyshev'
        Tens_alpha = alpha_T(A,B,C,arg);
    case 'Chebyshev2'
%         Tens_alpha = alpha_U(A,B,C);
        Tens_alpha = alpha_U_renorm(A,B,C,arg);
    case 'Taylor'
        Tens_alpha = alpha_Taylor(A,B,C,arg);
    case 'Gegenbauer'
%         Tens_alpha = alpha_C(A,B,C,arg);
        Tens_alpha = alpha_C_renorm(A,B,C,arg);
    otherwise
        error('This choice of basis is not implemented')
end

Tens = zeros(N,N,N);
if exist('intval','file') && isintval(arg)
    Tens = intval(Tens);
end
for m = 0:N-1
    for q = 0:N-1
        l = abs(q-m):2:min(q+m,N-1);
        idx = sub2ind([N,N,N],(q+1)*ones(1,length(l)),l+1,(q+l-m)/2+1);
        Tens(m+1,q+1,l+1) = Tens_alpha(idx);
    end
end

% SANITY CHECK
sum2one = max(max(sum(Tens,1)));
if sum2one > 1 + 1e-5
    warning("The computation of the linearization coefficients seems to be somewhat inaccurate")
    fprintf("\nmax_{m,n} sum_k alpha^k_{m,n} = %g, although it should be equal to 1\n",i2f(sum2one))
    pause
else
    save(str,'Tens')
    fprintf("\nLinearization coefficients computed and stored\n\n")
end

