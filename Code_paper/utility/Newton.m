function X = Newton(X,F,DF,para,it_max,tol,show)

sz = size(X);
    
res = tens2vect(F(X,para));
err = sum(abs(res));
if show
    display(err)
end
it = 0;
while (err > tol) && (it < it_max) && (err < 10^10)
    DX = vect2tens( DF(X,para) \ res, sz );
    X = X - DX;
    if isfield(para,'symmetrize')
        X = para.symmetrize(X);
    end
    res = tens2vect(F(X,para));
    err = sum(abs(res));
    if show
        display(err)
    end
    it = it+1;
end

if err > tol || isnan(err)
    warning('\nNewton method may not have converged, the residual error is %e',err)
end
