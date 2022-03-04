function plot_SH(u,L,nb,col,lw)

if nargin < 5
    lw = 1;
    if nargin < 4
        col = 'b';
        if nargin < 3
            nb = 1000;
            if nargin < 2
                L = pi;
            end
        end
    end
end

K = length(u);
pts = transpose(0:L/nb:L);

Mat = cos(pts*(0:K-1)*pi/L);
Mat(:,2:end) = 2*Mat(:,2:end);

plot(pts, Mat*u, 'color', col, 'Linewidth', lw)
xlabel('x')
ylabel('u')