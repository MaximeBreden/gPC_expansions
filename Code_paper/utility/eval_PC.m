function X_p = eval_PC(X, p, PC_type, mu)

[~,N] = size(X);

switch PC_type
    case 'Legendre'
        phi_p = legendreP(0:N-1,p);
    case 'Chebyshev'
        phi_p = chebyshevT(0:N-1,p);
    case 'Chebyshev2'
        phi_p = myChebyshevU(0:N-1,p); 
    case 'Taylor'
        phi_p = p.^(0:N-1);
    case 'Gegenbauer'
        phi_p = myGegenbauerC(0:N-1,mu,p);
    otherwise
        error('This choice of basis is not implemented')
end

X_p = X*transpose(phi_p);