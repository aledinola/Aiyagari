function [mu1,flag_mu,iter_final] = compute_mu_vec1(a_grid,pol_ap,pi_z,par)

% DESCRIPTION
% compute_mu_vec computes the stationary distribution mu(a,z) given policy
% function a'=pol_ap(a,z) and the Markov chain for exogenous shock
% pi_z(z,z'). Build big tran over (a,z). Uses a different method wrt
% compute_mu_vec.
% INPUTS
%   a_grid  : Fixed grid for assets a
%   pol_ap  : Policy function a'(a,z) where a' can be off grid
%   pi_z    : Transition matrix pi_z(z,z'), rows sum to one
%   par     : Struct with model parameters
% OUTPUTS
%   mu1        : Stationary distribution mu(a,z)
%   flag_mu    : Convergence flag, negative values mean no convergence
%   iter_final : Total no. of iterations
%
% Written by Alessandro Di Nola
% aledinola13@gmail.com

% Unpack
tol_dist     = par.tol_dist;
maxiter_dist = par.maxiter_dist;

% Set convergence flag to 0. If convergence fails, the flag will be changed
% to some negative value
%flag_mu = 0;

% STEP 1 - Find interpolating indexes and weights
% aInd is the location of the left point
[aInd,omega] = find_loc_vec(a_grid,pol_ap);

% STEP 2 - Call an external function to compute the distribution
[~,mu1,iter_final,flag_mu] = invDistrib(pi_z,aInd,omega,tol_dist,maxiter_dist);

sum_mu = sum(mu1,'all');
if abs(sum_mu-1)>1e-7
    error('MU does not sum to one!')
end

end %end function "compute_mu_vec1"

function [TEndog,invD,iter_final,conv_flag] = invDistrib(TExog,IndxOpt,pLowOpt,tol_dist,maxiter_dist)
% invDistrib computes the stationary distribution invD(endo,exo)
% INPUTS
%   TExog   :     nExog*nExog Markov chain exogenous shock
%   IndxOpt :     nEndog*nExog matrix with left points for a' (endo)
%   pLowOpt :     Weights for left points for a' (endo)
%   tol_dist:     Tolerance criterion
%   maxiter_dist: Maximum no. of iterations
% OUTPUTS
%   TEndog:   Large transition matrix
%   invD  :   Stationary distribution on (endo,exo)
% written by M.Reiter
% Modified by A. Di Nola on August 2022.
conv_flag  = 0;
iter_final = 0;
[nEndog,nExog] = size(IndxOpt);
nn = nEndog*nExog;
% For each exogenous state there is a different endogenous transition:
TEndog = sparse(nn,nn);
ee = (1:nEndog)';
for j=1:nExog
    ii = (j-1)*nEndog+1 : j*nEndog; % index of endogenous grid points
    % endogenous transition as it arises from linear interpolation
    jj   = IndxOpt(:,j);
    pLow = pLowOpt(:,j);
    TEndog(ii,ii) = sparse([jj;jj+1],[ee;ee],[pLow;1-pLow],nEndog,nEndog);
end
assert(all(abs(sum(TEndog)-1)<1e-10))

% start iteration with equal distribution:
D = ones(nn,1) / nn;
for iter = 1:maxiter_dist % maximum of 10000 iterations
    Dmid = TEndog*D;
    D2 = vec(reshape(Dmid,nEndog,nExog) * TExog);
    sumD = sum(D2);
    assert(abs(sumD-1)<1e-12);
    D2 = D2 / sumD;
    dist = max(abs(D-D2));
    D = D2;
    if(dist<tol_dist)
        invD = reshape(D,nEndog,nExog);
        iter_final = iter;
        return
    end
end
fprintf(1,'distance = %e',dist)
warning('maximum iterations exceeded')
invD = reshape(D,nEndog,nExog);
conv_flag = -1;

end % end function "invDistrib"

