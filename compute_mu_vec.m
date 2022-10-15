function [mu1,flag_mu,iter_final] = compute_mu_vec(a_grid,pol_ap,pi_z,par)

% DESCRIPTION
% compute_mu_vec computes the stationary distribution mu(a,z) given policy
% function a'=pol_ap(a,z) and the Markov chain for exogenous shock
% pi_z(z,z'). Build big tran over (a,z).
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
%   fun_genBigTran: Matlab function that computes the big transition matrix 
%                   on the joint Markov process (a,z). 
%
% Written by Alessandro Di Nola
% aledinola13@gmail.com

[na,nz] = size(pol_ap);

% Unpack
alg_ergo     = par.alg_ergo;
tol_dist     = par.tol_dist;
maxiter_dist = par.maxiter_dist;

% Initial condition for distribution: set in ergodicdist
%mu0  = ones(1,na*nz)/(na*nz); % this is arbitrary, as long as it sums to one
%mu0 = zeros(1,na*nz);
%mu0(1,1) = 1;

flag_mu = 0;

% Precompute indeces and weights required for interpolation of the assets
% policy function. IMPORTANT NOTE: the policy function a'(a,h,z) is fully
% summarized by aInd_arr(a,h,z) and omega_arr(a,h,z).

[Qmat] = fun_genBigTran(a_grid,pol_ap,pi_z);

[mu1,dist,iter_final] = ergodicdist(Qmat,tol_dist,maxiter_dist,alg_ergo);

% Find a better initial condition (otherwise we use mu0 defined above)
%[mu0,~]=eigs(Qmat',1);
%mu0=mu0'./sum(mu0);

% while dist>tol_dist && iter<=maxiter_dist
%     iter = iter+1;
%    
%     mu1 = mu0*Qmat;
%     
%     %Check
%     %check = sum(mu1(:));
%     
%     %Compute error
%     dist = max(abs(mu0-mu1));
%     % Update mu
%     mu0 = mu1;
%     
%     if disp_mu==1 && mod(iter,50)==0 
%         % Display some info every 50 iterations
%         %fprintf('sum(mu1) = %f\n', check);
%         fprintf('iter = %d, dist = %.10f\n', iter,dist);
%     end
%     
% end %end while
% iter_final = iter;

% Reshape to (na,nz) matrix
mu1 = reshape(mu1,na,nz);

sum_mu = sum(mu1,'all');
if abs(sum_mu-1)>1e-7
    error('MU does not sum to one!')
end

if dist>tol_dist
    flag_mu = -1;
end

end %end function "compute_mu_vec"

