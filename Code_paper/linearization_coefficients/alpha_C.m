function v=alpha_C(n,m,k,mu)

%Linearization coefficient \alpha^{(m,n)}_k for the Gegenbauer polynomials
%of parameter mu

Ind=(k<=min(m,n));
k=k.*Ind;
v=binbis(mu,k)./binbis(2*mu,n+m-2*k).*binbis(mu,n-k).*binbis(mu,m-k).*binbis(2*mu,m+n-k)./(binbis(mu,n+m-k)).*(n+m-2*k+mu)./(n+m-k+mu).*Ind;

%sqrt(binbis(2*lambda,n+m-2*k)./binbis(2*lambda,n)./binbis(2*lambda,m).*((n+lambda).*(m+lambda))./(lambda*(n+m-2*k+lambda)))...
%normalization factor for h_n=1 
    
function v=binbis(z,k)

%v=gamma(k+z)./(gamma(z).*factorial(k)); 
%The formulation below is better suited for large k and/or z, but does not 
%work for z<0.
v=ones(size(k));
Ind_k=k~=0;
v(Ind_k)=1./(k(Ind_k).*beta(z,k(Ind_k)));

end

end