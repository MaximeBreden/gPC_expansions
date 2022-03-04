function YY = resize(XX,Knew,Nnew)

[Kold,Nold] = size(XX);
Kold = (Kold+2)/6;

Kmin = min(Kold,Knew);
Nmin = min(Nold,Nnew);

Ind_K = (-Kmin+1:Kmin-1);
Ind_N = (1:Nmin);

YY = zeros(6*Knew-2,Nnew);
if  exist('intval','file') && isintval(XX(1))
    YY = intval(YY);
end
YY(1,Ind_N) = XX(1,Ind_N);
YY(Knew+1+Ind_K,Ind_N) = XX(Kold+1+Ind_K,Ind_N);
YY(3*Knew+Ind_K,Ind_N) = XX(3*Kold+Ind_K,Ind_N);
YY(5*Knew-1+Ind_K,Ind_N) = XX(5*Kold-1+Ind_K,Ind_N);
