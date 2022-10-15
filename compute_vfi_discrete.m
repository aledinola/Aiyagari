function [sol,flag,iter_final] = compute_vfi_discrete(par)

% DESCRIPTION
% compute_vfi computes the value function V(a,z) and policy functions for
% consumption and next-period assets c=pol_c(a,z), a'=pol_ap(a,z). They are
% all stored in the struct sol.
% INPUTS
%   par        : Struct with model parameters
% OUTPUTS
%   sol        : Struct with value and policy functions, all dim (na,nz)
%   flag       : Convergence flag, negative values mean no convergence
%   iter_final : Total no. of iterations
% AUXILIARY
%   sub_vfi_onestep : Updates the value function by max over a'.
%   sub_howard      : Howards policy improvement algorithm.
%
% Written by Alessandro Di Nola
% aledinola13@gmail.com

% Unpack
na     = par.na;
nz     = par.nz;
z_grid = par.z_grid;
a_grid = par.a_grid;
pi_z   = par.pi_z;
beta   = par.beta;
% tol_v  = par.tol_v;
% maxit_v = par.maxit_v;
% do_howard  = par.do_howard;
% disp_vfi   = par.disp_vfi;
% alg_howard = par.alg_howard;

opt(1) = par.tol_v; % Tolerance for VFI
opt(2) = par.maxit_v; % Maximum number of iterations (>1) for VFI
opt(3) = par.alg_lin;  % Algorithm (1: sparse LU, 2: BiConjugate Gradients
                          %                       Stabilized method)

[pol_ap_ind,V1,flag] = discdynprog('ReturnFn',a_grid,z_grid,pi_z,beta,opt,par);

pol_ap = a_grid(pol_ap_ind);

pol_c = zeros(na,nz);
for z_c = 1:nz
    a_val      = a_grid;
    z_val      = z_grid';
    aprime_val = pol_ap;
    [~,pol_c]  = ReturnFn(a_val,z_val,aprime_val,par);
end

if any(pol_c<0,'all')
    error('Some element of pol_c is negative!')
end

sol.V1     = V1;
sol.pol_ap = pol_ap;
sol.pol_c  = pol_c;

iter_final  = 1;

end %end function "compute_vfi_discrete"

