function [sol,flag,iter_final] = compute_vfi(par)

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
tol_v  = par.tol_v;
maxit_v = par.maxit_v;
do_howard  = par.do_howard;
disp_vfi   = par.disp_vfi;
alg_howard = par.alg_howard;

% Precompute cash on hand (RHS of the budget constraint)
cash_arr = zeros(na,nz);
for z_c = 1:nz
    z_val = z_grid(z_c);
    for a_c = 1:na
        a_val = a_grid(a_c);
        cash_arr(a_c,z_c) = fun.fun_cash(a_val,z_val,par);
    end
end


% Set initial condition for value function
V0 = zeros(na,nz);

flag = 1;

iter = 1;
err  = tol_v+1;

while err>tol_v && iter<=maxit_v
    
    [V1,pol_ap] = sub_vfi_onestepx(V0,cash_arr,par);
    
    if do_howard==1 && iter>2
        switch alg_howard
            case 1
                V1 = sub_howard(V1,pol_ap,pi_z,cash_arr,a_grid,par);
            case 2
                V1 = sub_howard_vec(V1,pol_ap,pi_z,cash_arr,a_grid,par);
            otherwise
                error('alg_howard is not specified correctly!')     
        end
    end
    
    % Check distance and update
    err = max(abs(V1(:)-V0(:)));
    V0  = V1;
    
    if disp_vfi==1
        fprintf('iter = %d, err = %f \n',iter,err)
    end
    iter = iter+1;
    
end
iter_final = iter;

if err>tol_v
    flag=-1;
end

% Compute other policy functions
pol_c = zeros(na,nz);
for z_c = 1:nz
    for a_c = 1:na
        cash = cash_arr(a_c,z_c);
        aprime = pol_ap(a_c,z_c);
        pol_c(a_c,z_c) = cash-aprime;
    end
end

if any(pol_c<0,'all')
    error('Some element of pol_c is negative!')
end

sol.V1     = V1;
sol.pol_ap = pol_ap;
sol.pol_c  = pol_c;

end %end function "compute_vfi"

