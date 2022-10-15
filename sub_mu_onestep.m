function [mu1] = sub_mu_onestep(mu,pi_z,aInd_arr,omega_arr)

% DESCRIPTION:
%   sub_mu_onestep does the distribution update.
% INPUTS:
%   "mu"         Initial distribution mu(a,z)
%   "pi_z"       Markov chain for exogenous shock z
%   "aInd_arr"   Array of interpolation indeces (a,z)
%   "omega_arr"  Array of interpolation weights (a,z)
%
% OUTPUTS:
%   "mu1"        Updated distribution mu1(a,z)
%
% NOTES:
% (1) For speed reasons, we have to precompute the interp indeces and
%     weights before calling this function (see input arguments aInd_arr,
%     omega_arr).
% (2) We used a speed-up trick suggested by Grey Gordon.
%     https://sites.google.com/site/greygordon/teaching
%     Minicourse on heterogeneous agent models
%-------------------------------------------------------------------------%

[na,nz] = size(mu); %State space: (a,z)

if ~isequal([na,nz],size(omega_arr))
    error('omega_arr is not compatible with mu')
end

% Initialize the "before-the-shocks" distribution mu_hat: after choice of a'
% is made but before the realization of z'.
mu_hat = zeros(na,nz);

for z_c = 1:nz
    for a_c = 1:na
        aInd       = aInd_arr(a_c,z_c);
        omega      = omega_arr(a_c,z_c);
        % NOTE: at this stage, we don't multiply by pi(z,z')
        mu_hat(aInd,z_c)   = mu_hat(aInd,z_c)+omega*mu(a_c,z_c);
        mu_hat(aInd+1,z_c) = mu_hat(aInd+1,z_c)+(1-omega)*mu(a_c,z_c);
    end %a
    
end %z

% Matrix multiplication mu_hat(a',z)*pi_z(z,z')==> mu(a',z')
mu1 = mu_hat*pi_z;

end %end function "sub_mu_onestep"

