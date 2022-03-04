function plot_Lorenz(u,nb,col,f1,f2)

K = (length(u)+2)/6;

Omega = real(u(1));
x = u(2:2*K);
y = u(2*K+1:4*K-1);
z = u(4*K:6*K-2);

T = 2*pi/Omega;
t = transpose(0:T/nb:T);
M = exp(1i*Omega*t*(-K+1:K-1));

valx = real(M*x);
valy = real(M*y);
valz = real(M*z);

if nargin < 3
    col = 'b';
    if nargin < 2
        nb = 500;
    end
end

if nargin < 4
    figure(1)
else
    figure(f1)
end
plot3(valx,valy,valz,col,'linewidth',2)
xlabel('x')
ylabel('y')
zlabel('z')
hold on

if nargin < 5
    figure(2)
else
    figure(f2)
end
T_resc = 2*pi;
t_resc = transpose(0:T_resc/nb:T_resc);

s1_resc = subplot(3,1,1);
plot(t_resc,valx,col)
xlabel('t (rescaled)')
ylabel('x(t)')
axis tight
hold(s1_resc,'on')

s2_resc = subplot(3,1,2);
plot(t_resc,valy,col)
xlabel('t (rescaled)')
ylabel('y(t)')
axis tight
hold(s2_resc,'on')

s3_resc = subplot(3,1,3);
plot(t_resc,valz,col)
xlabel('t (rescaled)')
ylabel('z(t)')
axis tight
hold(s3_resc,'on')
