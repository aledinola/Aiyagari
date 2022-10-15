function [Qmat] = fun_genBigTran(a_grid,pol_ap,pi_z)

% DESCRIPTION:
%   fun_genBigTran generates the big Markov chain "Qmat" over the entire
%   state-space (a,z).
% INPUTS:
%   "a_grid"   Exogennous grid for assets, dim: (na,1)
%   "pol_ap"   Policy function a'(a,z) where a' can be off grid
%   "pi_z"     Markov chain for exogenous shock z
%
% OUTPUTS:
%   "Qmat"     Big transition matrix on (a,z), dim: (na*nz,na*nz)
%
% AUXILIARY:
%   "find_loc_vec"  Used to find interp indeces and weights (based on histc)
%
% Written by Alessandro Di Nola
% aledinola13@gmail.com

[na,nz] = size(pol_ap);
NA = (1:na)';

% STEP 1 - Find interpolating indexes and weights
% aInd is the location of the left point
[aInd,omega] = find_loc_vec(a_grid,pol_ap);

% STEP 2 - For each z, build an na*na matrix transition matrix, called G_z, 
% that represents the policy function a'=g(a,z). Version with loops is 
% commented below.
G = cell(nz,1);
for z_c = 1:nz
   G{z_c} = sparse(NA,aInd(:,z_c),omega(:,z_c),na,na)+...
       sparse(NA,aInd(:,z_c)+1,1-omega(:,z_c),na,na);
end
% G is (nz*na,na) matrix

% G = cell(nz,1);
% for z_c = 1:nz
%    G{z_c} = sparse(na,na);
%    for a_c = 1:na
%        [aInd,omega] = find_loc(a_grid,pol_ap(a_c,z_c));
%        G{z_c}(a_c,aInd)   = omega;
%        G{z_c}(a_c,aInd+1) = 1-omega;
%    end
% end
% % G is (nz*na,na) matrix

% STEP 3 - Introduce the transition probs for the exogenous shock z
Q = cell(nz,1);
for z_c = 1:nz
    % Kron product of (1,nz) with (na,na) gives (na,nz*na) matrix
    Q{z_c} = kron(pi_z(z_c,:),G{z_c}); %dim: (na,na*nz)
end %close z
% Note: vertcat with {:} operator seems faster than cell2mat.
%Qmat = cell2mat(Q); %dim: (na*nz,na*nz), rows sum to one
Qmat = vertcat(Q{:});

end %end function "fun_genBigTran"

