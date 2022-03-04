function XX = symmetrize(XX)

XX_old = XX;
[K,~] = size(XX);
K = (K+2)/6;

x = XX(2:2*K,:);
y = XX(2*K+1:4*K-1,:);
z = XX(4*K:6*K-2,:);

XX(1,:) = real(XX(1,:));
xsym = ( x(K:2*K-1,:) +  conj(flip(x(1:K,:))) ) / 2;
XX(2:2*K,:) = [conj(flip(xsym(2:K,:))); xsym];
ysym = ( y(K:2*K-1,:) +  conj(flip(y(1:K,:))) ) / 2;
XX(2*K+1:4*K-1,:) = [conj(flip(ysym(2:K,:))); ysym];
zsym = ( z(K:2*K-1,:) +  conj(flip(z(1:K,:))) ) / 2;
XX(4*K:6*K-2,:) = [conj(flip(zsym(2:K,:))); zsym];

if sum(sum(abs(XX-XX_old)))>10^-6
    warning("The numerical approximation seems far from real")
end