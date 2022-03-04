function YY = resize_deter(XX,Knew)

Kold = (length(XX)+2)/6;
Kmin = min(Kold,Knew);
Ind = (-Kmin+1:Kmin-1);
YY = zeros(6*Knew-2,1);
YY(1) = XX(1);
YY(Knew+1+Ind) = XX(Kold+1+Ind);
YY(3*Knew+Ind) = XX(3*Kold+Ind);
YY(5*Knew-1+Ind) = XX(5*Kold-1+Ind);
