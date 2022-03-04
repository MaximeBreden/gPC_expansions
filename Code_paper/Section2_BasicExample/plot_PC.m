function plot_PC(x, PC_type, mu, col, nb_pts, f, lw)

if nargin < 7
    lw = 1;
    if nargin <6
        figure
        if nargin < 5
            nb_pts = 100;
            if nargin < 4
                col = 'b';
                if nargin < 3
                    mu = 0;
                end
            end
        end  
    else
        figure(f)
    end
end

if size(x,2) == 1
    x = transpose(x);
end

N = length(x);

if N <= 40 %this is the easiest way to plot the approximate solution, but seems to be unstable for N > 40
    syms p
    fplot(eval_PC(x, p, PC_type, mu), [-1,1], col, 'Linewidth', lw)
else %stabler but slower
    fprintf("\nPlotting, this may take a couple of seconds ...\n")
    p = linspace(-1,1,nb_pts);
    phi_p = zeros(1,nb_pts);
    for i = 1:nb_pts
        phi_p(i) = eval_PC(x, p(i), PC_type, mu);
    end
    plot(p, phi_p, col, 'Linewidth', lw)
end