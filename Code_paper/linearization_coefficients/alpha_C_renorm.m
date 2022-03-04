function alpha = alpha_C_renorm(n,m,k,mu)

%Linearization coefficient \alpha^{(m,n)}_k for the Gegenbauer polynomials
%of parameter mu

Ind = (k<=min(m,n));
k = k.*Ind;
if nargin > 3 && exist('intval','file') && isintval(mu)
   k = intval(k);
   m = intval(m);
   n = intval(n);
end
alpha = binbis(mu,n-k) ./ binbis(2*mu,n) .* binbis(mu,m-k) ./ binbis(2*mu,m) .* binbis(mu,k) .* binbis(2*mu,m+n-k) ./ (binbis(mu,n+m-k)) .* (n+m-2*k+mu) ./ (n+m-k+mu) .* Ind;
    
function v = binbis(z,k)
%v=gamma(k+z)./(gamma(z).*factorial(k)); 
%The formulation below is better suited for large k and/or z, but does not 
%work for z<0.
v = ones(size(k));
if exist('intval','file') && isintval(k(1))
    v = intval(v);
end
Ind_k = (k~=0);
v(Ind_k) = 1./(k(Ind_k) .* beta(z,k(Ind_k)));
end

end