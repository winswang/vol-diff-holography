function x = normMax(x, mode)
if nargin < 2
    mode = 0;
end

xmax = max(x(:));
if mode == 0
    x = x/xmax;
else
    xmin = min(x(:));
    x = (x - xmin)/(xmax - xmin);
end