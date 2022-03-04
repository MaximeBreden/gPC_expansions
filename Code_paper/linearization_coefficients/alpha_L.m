function alpha = alpha_L(n,m,k,arg)

%Linearization coefficient \alpha^{(m,n)}_k for the Legendre polynomials

% alpha=0*n;
% for i=1:numel(k)
%     if k(i)>min(n(i),m(i))
%         alpha(i)=0;
%     else
%         alpha(i)=binbis(1/2,k(i)).*binbis(1/2,n(i)-k(i)).*binbis(1/2,m(i)-k(i))./binbis(1/2,n(i)+m(i)-k(i)).*(n(i)+m(i)-2*k(i)+1/2)/(n(i)+m(i)-k(i)+1/2);
%     end
% end

Ind = (k<=min(m,n));
k=k.*Ind;
if nargin > 3 && exist('intval','file') && isintval(arg)
   k = intval(k);
   m = intval(m);
   n = intval(n);
end
alpha = binbis(1/2,k) .* binbis(1/2,n-k) .* binbis(1/2,m-k) ./ binbis(1/2,n+m-k) .* (n+m-2*k+1/2) ./ (n+m-k+1/2) .* Ind;


function v = binbis(z,k)
%v=gamma(k+z)./(gamma(z).*factorial(k));
v = ones(size(k));
if exist('intval','file') && isintval(k(1))
    v = intval(v);
end
Ind_k = (k~=0);
v(Ind_k) = 1./(k(Ind_k) .* beta(z,k(Ind_k)));
end

end