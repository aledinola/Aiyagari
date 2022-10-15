function [V] = sub_howard_vec(~,pol_ap,pi_z,cash_arr,a_grid,par)

% DESCRIPTION
% sub_howard_vec implements the Howard's policy improvement algorithm as
% outlined in LS (add ref). This version relies on calculating the large
% transition matrix over the joint process (a,z). While faster than
% sub_howard, it is much more memory-intensive.
% INPUTS
%   V          : Value function (OPTIONAL, not used)
%   pol_ap     : Policy function a'(a,z) where a' can be off grid
%   pi_z       : Transition matrix pi_z(z,z'), rows sum to one
%   cash_arr   : Precomputed cash-on-hand, dim: (na,nz)
%   a_grid     : Assets grid, dim: (na,1)
%   par        : Struct with model parameters
% OUTPUTS
%   V          : Value function (inout)
% AUXILIARY
%   fun_genBigTran: Matlab function that computes the big transition matrix 
%                   on the joint Markov process (a,z), dim:(na*nz,na*nz).
% Written by Alessandro Di Nola
% aledinola13@gmail.com

% Unpack
beta    = par.beta;
gamma   = par.gamma;
alg_lin = par.alg_lin;
[na,nz] = size(pol_ap);

% STEP 1 - Precompute static payoff at optimal policies
cons   = cash_arr-pol_ap;  % (na,nz)
payoff = fun.util(cons,gamma); % (na,nz)
payoff_vec = reshape(payoff,na*nz,1);

% STEP 2 - Compute big transition matrix Q(x,x') where x = (a,z)
% a varies first, so a(1),z(1),a(2),z(1),..,a(na),z(1),a(1),z(2), etc.
Qmat = fun_genBigTran(a_grid,pol_ap,pi_z);

% STEP 3 - Solve large linear system
% V = payoff + beta*Qmat*V ==> (I-beta*Qmat)*V = payoff
% This can be solved as A*x=b where A = I-beta*Qmat and b=payoff
% where
% V is (na*nz,1),
% payoff is (na*nz,1),
% beta is a scalar,
% Qmat us (na*nz,na*nz)

Ibig = speye(na*nz,na*nz);
switch alg_lin
    
    case 1
        V_vec = (Ibig-beta*Qmat)\payoff_vec;
    case 2
        V_vec = sparse(na*nz,1); %zeros(na*nz,1);
        metric = 10;
        while metric>1e-8
           V_vec1 = payoff_vec + beta*Qmat*V_vec;
           metric = max(abs(V_vec1-V_vec));
           V_vec  = V_vec1;
        end
        
        %[V_vec,~] = bicgstab((Ibig-beta*Qmat),payoff_vec);
        %if flag>0
        %    warning('bicgstab did not converge')
        %end
    otherwise
        error('alg_lin not specified correctly!')
        
end

V = reshape(V_vec,na,nz);

end %end sub_howard_vec

