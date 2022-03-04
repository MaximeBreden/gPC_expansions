function val = myGegenbauerC(n,mu,x)

val = gegenbauerC(n,mu,x)./binbis(2*mu,n);

function v=binbis(z,k)

%v=gamma(k+z)./(gamma(z).*factorial(k)); 
%The formulation below is better suited for large k and/or z, but does not 
%work for z<0.
v=ones(size(k));
Ind_k=k~=0;
v(Ind_k)=1./(k(Ind_k).*beta(z,k(Ind_k)));

end

end