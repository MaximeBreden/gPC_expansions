function DF = DF_SH_deter_cont(u, para)

% Derivative of the map F with restpect to (rho,u)

DF = [u, DF_SH_deter(u, para)];