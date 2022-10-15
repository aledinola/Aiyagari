function [jl,omega] = find_loc(x_grid,xi)

% DESCRIPTION:
% Find jl s.t. x_grid(jl)<=xi<x_grid(jl+1)
% for jl=1,..,N-1
% omega is the weight on x_grid(jl)
% ASSUMPTIONS:
% x_grid must be a column vector
% xi must be scalar
% If these assumptions are not met, the routine will not work as intended.
% NOTES:
% This function is not vectorized; this means that input xi
% must be a scalar. For a vectorized version, see find_loc_vec
% or find_loc_vec1.

nx = size(x_grid,1);

jl = max(min(locate(x_grid,xi),nx-1),1);
%Weight on x_grid(j)
omega = (x_grid(jl+1)-xi)/(x_grid(jl+1)-x_grid(jl));
omega = max(min(omega,1),0);

end %end function "find_loc"