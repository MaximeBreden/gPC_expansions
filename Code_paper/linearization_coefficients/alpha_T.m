function alpha = alpha_T(n,m,k,arg)

%Linearization coefficient \alpha^{(m,n)}_k for the Chebyshev polynomials
%of the first kind

% alpha=0*n;
% for i=1:numel(n)
%     if k(i)==0 || k(i)==min(n(i),m(i))
%         if min(n(i),m(i))==0
%             alpha(i)=1;
%         else
%             alpha(i)=1/2;
%         end
%     end
% end

k0 = (k==0);
Min = min(m,n);
Min0 = (Min==0);

alpha = k0.*Min0 + 1/2*(1-Min0).*max(k0,k==Min);
