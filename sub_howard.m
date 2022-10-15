function [V] = sub_howard(V,pol_ap,pi_z,cash_arr,a_grid,par)

% DESCRIPTION
% sub_howard implements the Howard's policy improvement algorithm as
% outlined in LS (add ref). This version is mainly loop-based: while slower
% than sub_howard_vec, it is less memory-intensive. Moreover, if the
% state-space is large (i.e. if na*nz is a large number), it might be even
% faster than sub_howard_vec.
% INPUTS
%   V          : Value function (inout)
%   pol_ap     : Policy function a'(a,z) where a' can be off grid
%   pi_z       : Transition matrix pi_z(z,z'), rows sum to one
%   cash_arr   : Precomputed cash-on-hand, dim: (na,nz)
%   a_grid     : Assets grid, dim: (na,1)
%   par        : Struct with model parameters
% OUTPUTS
%   V          : Value function (inout)
% AUXILIARY
%   interp1qr  : Matlab function to do interp1(X,Y,XI).
%   myinterp1q : MEX function to do interp1(X,Y,XI), written in Fortran.
% NOTES
%   It is a bit faster to use myinterp1q rather than interp1qr. Take into
%   account, however, that MEX-files are less stable than normal Matlab
%   functions.
% Written by Alessandro Di Nola
% aledinola13@gmail.com

% Unpack
beta  = par.beta;
gamma = par.gamma;
n_how = par.n_how;
[na,nz] = size(pol_ap);

% STEP 1 - Precompute static payoff at optimal policies. Version with loops
% is commented below.
cons   = cash_arr-pol_ap;  % (na,nz)
payoff = fun.util(cons,gamma); % (na,nz)

% payoff = zeros(na,nz);
% for z_c = 1:nz
%     for a_c = 1:na
%           cons = cash_arr(a_c,z_c)-pol_ap(a_c,z_c);
%           payoff(a_c,z_c) = util(cons,gamma);
%     end
% end

% STEP 2 - Find interpolating indexes and weights
% aInd is the location of the left point, dim: (na,nz)
% aInd is the weight the left point,      dim: (na,nz)
[aInd,omega] = find_loc_vec(a_grid,pol_ap);

% STEP 3 - Run Howards iterations. A simpler version is commented below.
V_update = zeros(na,nz);
for howard = 1:n_how
    EV = V*pi_z';
    for z_c = 1:nz
        EVz = EV(:,z_c);
        for a_c = 1:na
            V_update(a_c,z_c) = payoff(a_c,z_c)+beta*( omega(a_c,z_c)*EVz(aInd(a_c,z_c))+(1-omega(a_c,z_c)).*EVz(aInd(a_c,z_c)+1) );
        end
    end
    V = V_update;
end

% STEP 2 - Run Howards iterations. Version with loops is commented below.
% V_update = zeros(na,nz);
% for howard = 1:100
%     EV = V*pi_z';
%     for z_c = 1:nz
%         EVz = EV(:,z_c);
%         V_update(:,z_c) = payoff(:,z_c)+beta*myinterp1q(a_grid,EVz,pol_ap(:,z_c));
%         %V_update(:,z_c) = payoff(:,z_c)+beta*( omega(:,z_c).*EVz(aInd(:,z_c))+(1-omega(:,z_c)).*EVz(aInd(:,z_c)+1) );
%     end
%     V = V_update;
% end

end %end sub_howard

