function y = my_beta(z,w)

%Usable even for negative inputs, which is not the case for Matlab built-in
%beta function.

y=gamma(z)*gamma(w)/gamma(z+w);