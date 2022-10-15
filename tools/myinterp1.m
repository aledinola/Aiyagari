function yi = myinterp1(x,y,xi)

% PURPOSE:
% myinterp1 performs 1-D linear interpolation (with extrapolation). For
% speed reasons it does not do any input checking.
%
% INPUTS:
%   "x"    Grid, must be a sorted column vector, dim: (nx,1)
%   "y"    Column vector with function values, dim: size(x)
%   "xi"   query point and MUST be a scalar
%   If xi<x(1) or xi>x(end), we do extrapolation.
% OUTPUT
%   "yi" is the value of the interpolated function at query point xi 
%
% -------------------------------------------------------------------------
% Copyright © 2018 by Alessandro Di Nola. All rights 
% reserved.
% -------------------------------------------------------------------------

n = size(x,1);
if size(y,1)~=n
    error('myinterp1: x and y must be vectors with same length!')
end

% STEP 1 - Find x(j)<= xi <x(j+1), for j=1,..,n-1
j = locate(x,xi);
% xi is between x(j) and xx(j+1)
%
% j = 0 and j = n means x is out of range!!!!!!!!!!!!!!

j = max(min(j,n-1),1);

% STEP 2 - Apply linear interpolation formula
slope = (y(j+1)-y(j))/(x(j+1)-x(j));
yi = y(j)+(xi-x(j))*slope;

% % Give NaN to the values of 'yi' corresponding to 'xi' out of the range of 'x'
% if extrap==0
%     yi(xi<x(1) | xi>x(end)) = NaN;
% end

end %end function