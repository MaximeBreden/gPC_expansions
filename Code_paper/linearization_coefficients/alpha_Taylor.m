function alpha = alpha_Taylor(n,m,k,arg)

%Linearization coefficient \alpha^{(m,n)}_k for the monomial polynomials

% alpha=0*n;
% for i=1:numel(n)
%     if k(i)==0
%         v(i)=1;
%     end
% end

alpha = double(k==0);
