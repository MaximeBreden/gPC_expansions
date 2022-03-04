function X_p = eval_PC_2para(X, p, PC_type, mu)

[K,N1,N2] = size(X);

switch PC_type{1}
    case 'Legendre'
        phi_p1 = legendreP(0:N1-1,p(1));
    case 'Chebyshev'
        phi_p1 = chebyshevT(0:N1-1,p(1));
    case 'Chebyshev2'
        phi_p1 = myChebyshevU(0:N1-1,p(1)); 
    case 'Taylor'
        phi_p1 = p(1).^(0:N1-1);
    case 'Gegenbauer'
        phi_p1 = myGegenbauerC(0:N1-1,mu(1),p(1));
    otherwise
        error('This choice of basis is not implemented')
end

switch PC_type{2}
    case 'Legendre'
        phi_p2 = legendreP(0:N2-1,p(2));
    case 'Chebyshev'
        phi_p2 = chebyshevT(0:N2-1,p(2));
    case 'Chebyshev2'
        phi_p2 = myChebyshevU(0:N2-1,p(2)); 
    case 'Taylor'
        phi_p2 = p(2).^(0:N2-1);
    case 'Gegenbauer'
        phi_p2 = myGegenbauerC(0:N2-1,mu(2),p(2));
    otherwise
        error('This choice of basis is not implemented')
end

phi_p = permute( repmat( transpose(phi_p1) * phi_p2, 1, 1, K), [3 1 2]);

X_p = sum(X.*phi_p, [2 3]);