function [jl,omega] = find_loc_vec(x_grid,xi)

% DESCRIPTION
% Inputs are column vector x_grid (nx,1) and the column vector xi (nxi,1).
% For each i=1,..,nxi, find jl(i) s.t. x_grid(jl(i))<=xi(i)<x_grid(jl(i)+1)
% for jl=1,..,nx-1. 
% Outputs:
% jl is integer (nxi,1) and omega is a double (nxi,1).
% omega is the weight on x_grid(jl)
% IMPORTANT
% x_grid must be a column vector, NOT a row vector.
% NOTES
% Similar to find_loc_vec1 but here we use histc instead of interp1. If the
% input xi is a scalar, find_loc is faster.
% Matlab does not recommend histc. Instead, we can use discretize
% (histcounts is too slow!!).
% It seems that this works even if xi is a matrix

%nx = size(x_grid,1);

% For each 'xi', get the position of the 'x' element bounding it on the left [p x 1]
[~,jl] = histc(xi,x_grid);
jl(xi<=x_grid(1))=1;
jl(xi>=x_grid(end))=length(x_grid)-1;


%Weight on x_grid(j)
omega = (x_grid(jl+1)-xi)./(x_grid(jl+1)-x_grid(jl));
omega = max(min(omega,1),0);

%validateattributes(omega, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0,'<=', 1})

end %end function "find_loc_vec"