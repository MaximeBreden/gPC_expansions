function alpha = alpha_U_renorm(n,m,k,arg)

%Linearization coefficient \alpha^{(m,n)}_k for the Chebyshev polynomials
%of the second kind, normalized so that U_n(1) = 1.

ind = double(k<=min(n,m));
if nargin > 3 && exist('intval','file') && isintval(arg)
   k = intval(k);
   m = intval(m);
   n = intval(n);
end
renorm = (m+n-2*k+1) ./((m+1).*(n+1));
alpha = ind .* renorm;