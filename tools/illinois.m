function x = illinois(f,x0,x1,xtol,ftol,maxiter)

%
% Purpose: Solve a single non-linear equation using the Illinois method
%
% Written by Paul Klein in August 2019
%
% f is a string containing the name of the function
% x0 is a lower bound of the solution
% x1 is an upper bound of the solution
% ftol is the function tolerance
% xtol is the solution tolerance
% maxiter is the maximum number of iterations allowed
% Modified by Alessandro Di Nola in August 2022

xdev = 1;
iter = 0;

y0 = feval(f,x0);
y1 = feval(f,x1);

x = x0;
y = y0;

while (xdev>xtol || abs(y)>ftol) && iter<maxiter
    newx = (x0*y1 - x1*y0)/(y1-y0);
    xdev = abs(x-newx);
    x = newx;
    y = feval(f,x);
    if y*y1<0
        x0 = x1;
        y0 = y1;
        x1 = x;
        y1 = y;
    else
        x1 = x;
        y1 = y;
        y0 = 0.5*y0;
    end
    iter = iter + 1;
end

end %end function
