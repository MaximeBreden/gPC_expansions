function y = i2f(x,opt)

%Routine for dealing with data without having to worry whether the input is
%a usual float or an intval. If the input is an intval, the default output
%is the supremum

if exist('intval','file') && isintval(x(1))
    if nargin == 2 
        switch opt
            case 'inf'
                y = inf(x);
            case 'sup'
                y = sup(x);
            case 'mid'
                y = mid(x);
        end
    else
        y = sup(x);
    end
else
    y = x;
end