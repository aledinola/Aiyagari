function [mu1,flag_mu,iter_final] = compute_mu(a_grid,pol_ap,pi_z,par)

% DESCRIPTION
%   compute_mu computes the stationary distribution mu(a,z) given policy
%   function a'=pol_ap(a,z) and the Markov chain for exogenous shock
%   pi_z(z,z'). It uses a loop-based method that, while slower, is less 
%   memory-intensive.
% INPUTS
%   a_grid  : Fixed grid for assets a
%   pol_ap  : Policy function a'(a,z) where a' can be off grid
%   pi_z    : Transition matrix pi_z(z,z'), rows sum to one
%   par     : Struct with model parameters
% OUTPUTS
%   mu1        : Stationary distribution mu(a,z)
%   flag_mu    : Convergence flag, negative values mean no convergence
%   iter_final : Total no. of iterations
% AUXILIARY
%   sub_mu_onestep: This function does one step of the distrib. update
%   find_loc      : Find interpolation indeces and weights
%
% Written by Alessandro Di Nola
% aledinola13@gmail.com


[na,nz] = size(pol_ap);

% Unpack
tol_dist = par.tol_dist;
maxiter_dist = par.maxiter_dist;
disp_mu = par.disp_mu;

% Initial condition for distribution
mu0  = ones(na,nz)/(na*nz); % this is arbitrary, as long as it sums to one

dist    = 10;
iter    = 0;
flag_mu = 0;

% Precompute indeces and weights required for interpolation of the assets
% policy function. IMPORTANT NOTE: the policy function a'(a,z) is fully
% summarized by aInd_arr(a,z) and omega_arr(a,z).
aInd_arr  = ones(na,nz);
omega_arr = zeros(na,nz);
for z_c = 1:nz    
    for a_c = 1:na
        aprime     = pol_ap(a_c,z_c);
        [aInd,omega] = find_loc(a_grid,aprime);
        aInd_arr(a_c,z_c)  = aInd;
        omega_arr(a_c,z_c) = omega;
    end
end

while dist>tol_dist && iter<=maxiter_dist
    iter = iter+1;
   
    mu1 = sub_mu_onestep(mu0,pi_z,aInd_arr,omega_arr);
    
    %Check
    check = sum(mu1(:));
    
    %Compute error
    dist = max(abs(mu0(:)-mu1(:)));
    % Update mu
    mu0 = mu1;
    
    if disp_mu==1 && mod(iter,50)==0 
        % Display some info every 50 iterations
        fprintf('sum(mu1) = %f\n', check);
        fprintf('iter = %d, dist = %.10f\n', iter,dist);
    end
    
end %end while
iter_final = iter;

sum_mu = sum(mu1,'all');
if abs(sum_mu-1)>1e-7
    error('MU does not sum to one!')
end

if dist>tol_dist
    flag_mu = -1;
else
    fprintf('MU (distr.) converged after %d iterations! \n', iter);
end

end %end function "compute_mu"

