function [tab_sol, tab_bif] = continuation_SH_pal(F, DF, pt0, para0, vect0, ds, depth, nb_pts_max, rho_max, tol, it_max)  

fprintf("\nCurrent depth: %i\n", depth)

tab_sol = pt0; % table where we store the equilibria (+ corresponding parameter value)
tab_bif = []; % table where we store the equilibria (+ corresponding parameter value) at which a pitchfork seems to occur

if isempty(vect0)
    % if there is no initial tangent vector given, we compute one
    vect0 = null(DF(u0,para0));
    if size(vect0,2) > 1
        warning("\nThe starting point seems to be singular\n")
    end
    vect0 = vect0(:,1);
end

DF0 = DF(pt0(2:end),para0);
det0 = det([vect0';DF0]);
nb_pts = 0;
while nb_pts < nb_pts_max && para0.rho < rho_max
    pt1_pred = pt0 + ds * vect0; % predictor
    para_pred = para0;
    para_pred.rho = pt1_pred(1);
    [pt1, para1] = Newton_pal(F, DF, pt1_pred, para_pred, pt1_pred, vect0, tol, it_max); % corrector using Newton's method
    DF1 = DF(pt1(2:end),para1);
    V = null(DF1); % computation of a new tangent vector at the new point pt1
    if size(V,2) == 1
        vect1 = V;
    else % close enough to a bifurcation that numerically the kernel becomes more than 1d
        fit = -1;
        for i = 1:size(V,2) % we try to select the "right" vector in the kernel to keep going in the same direction until we really cross the bifurcation
            vect_temp = V(:,i);
            psc = abs(vect_temp'*vect0);
            if psc > fit
               vect1 = vect_temp;
               fit = psc;
            end
        end
    end
    if vect1'*vect0 < 0 % checking that we still go in the same direction
        vect1 = - vect1;
    end
    det1 = det([vect1'; DF1]);
    if det0*det1 > 0 % no bifurcation, we store the newly computed point and keep going
        tab_sol = [tab_sol, pt1];
        nb_pts = nb_pts + 1;
        pt0 = pt1;
        para0 = para1;
        vect0 = vect1;
        DF0 = DF1;
        det0 = det1;
    else % bifurcation
        fprintf("\nBifurcation detected!\n")
        tab_bif = [tab_bif, (pt0+pt1)/2];
        if depth > 0 % we compute the bifurcating direction and restart on each branch
            [~,~,V] = svd(DF0); % the last columns of V should approximately span the kernel of the Jacobian DF0
            vect_cont = V(:,end); 
            vect_bif = V(:,end-1); 
            if abs(vect_cont'*vect0) < abs(vect_bif'*vect0) % checking that we correctly identified which vector is continuing the original branch and which one is bifurcating
                vect_temp = vect_cont;
                vect_cont = vect_bif;
                vect_bif = vect_temp;
            end
            if vect_cont'*vect0 < 0 % checking that the continuation vector still goes in the same direction
                vect_cont = -vect_cont;
            end
            
            pt10_pred = pt0 + ds * vect_cont; % predictor on the original branch
            para_pred10 = para0;
            para_pred10.rho = pt10_pred(1);
            [pt10, para10] = Newton_pal(F, DF, pt10_pred, para_pred10, pt10_pred, vect_cont, tol, it_max); % corrector on the original branch
            [tab_sol_10, tab_bif_10] = continuation_SH_pal(F, DF, pt10, para10, vect_cont, ds, depth-1, nb_pts_max, rho_max, tol, it_max); % continuation on the original branch

            pt11_pred = pt0 + ds * vect_bif; % predictor on one side of the bifurcating branch
            para_pred11 = para0;
            para_pred11.rho = pt11_pred(1);
            [pt11, para11] = Newton_pal(F, DF, pt11_pred, para_pred11, pt11_pred, vect_bif, tol, it_max); % corrector on that side of the bifurcating branch
            [tab_sol_11, tab_bif_11] = continuation_SH_pal(F, DF, pt11, para11, vect_bif, ds, depth-1, nb_pts_max, rho_max, tol, it_max); % continuation on that side of the bifurcating branch

            pt12_pred = pt0 - ds * vect_bif; % predictor on the other side of the bifurcating branch
            para_pred12 = para0;
            para_pred12.rho = pt12_pred(1);
            [pt12, para12] = Newton_pal(F, DF, pt12_pred, para_pred12, pt12_pred, -vect_bif, tol, it_max); % corrector on that other side of the bifurcating branch
            [tab_sol_12, tab_bif_12] = continuation_SH_pal(F, DF, pt12, para12, -vect_bif, ds, depth-1, nb_pts_max, rho_max, tol, it_max); % continuation on that other side of the bifurcating branch

            tab_sol = [tab_sol, tab_sol_10, tab_sol_11, tab_sol_12];
            tab_bif = [tab_bif, tab_bif_10, tab_bif_11, tab_bif_12];
            return
        else
            fprintf("\nMaximal depth reached on the current branch\n")
            return
        end
    end
end
if nb_pts == nb_pts_max 
    fprintf("\nMaximal number of points reached on the current branch\n")
else
    fprintf("\nMaximal value of rho reached on the current branch\n")
end

end

function [pt, para] = Newton_pal(F, DF, pt, para, p0, v0, tol, it_max)
%%% Newton method for the enlarged map containing F + the pseudo arclength continuation equation
res = [v0'*(pt-p0); F(pt(2:end),para)];
err = sum(abs(res));
it = 0;
while (err > tol) && (it < it_max)
    pt = pt - [v0';DF(pt(2:end),para)] \ res;
    para.rho = pt(1);
    res = [v0'*(pt-p0); F(pt(2:end),para)];
    err = sum(abs(res));
    it = it + 1;
end
if err > tol
    warning('\nNewton method may not have converged, the residual error is %e',err)
end
end
