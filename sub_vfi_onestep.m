function [V1,pol_ap] = sub_vfi_onestep(V0,cash_arr,par)

% DESCRIPTION
%   sub_vfi_onestep implements the Bellman operator: given an initial V0 it
%   performs the maximization and returns an updated V1. It uses the
%   Parallel computing toolbox.
% INPUTS
%   V0         : Value function (initial condition), dim: (na,nz)
%   cash_arr   : Precomputed cash-on-hand, dim: (na,nz)
%   par        : Struct with model parameters
% OUTPUTS
%   V1         : Updated value function, dim: (na,nz)
% AUXILIARY
%   golden     : this function does point-wise golden search method.
%
% Written by Alessandro Di Nola
% aledinola13@gmail.com

% Unpack
beta   = par.beta;
gamma  = par.gamma;
na     = par.na;
nz     = par.nz;
pi_z   = par.pi_z;
a_grid = par.a_grid;
tol_golden = par.tol_golden;

% Pre-allocate memory
V1     = zeros(na,nz);
pol_ap = zeros(na,nz);

% Set lower bound for a'
ap_lb = a_grid(1);
% Set upper bound for a'. This is a loose upper bound 
ap_ub = a_grid(na); 

% Compute expected value
EV = V0*pi_z'; % EV(a',z)

parfor z_c = 1:nz
    EVz = EV(:,z_c);
    for a_c = 1:na
        
        % Set lower and upper bound for a'. The upper bound cannot be
        % larger than max(a_grid); so we avoid doing extrapolation (risky!)
        cash  = cash_arr(a_c,z_c);
        ap_ub1 = min(ap_ub,cash-1e-10);
        
        fun_rhs = @(ap) fun.util(cash-ap,gamma) + beta*myinterp1(a_grid,EVz,ap);
        
        % Maximize the RHS of the Bellman equation
        optimal_aprime  = golden(fun_rhs,ap_lb,ap_ub1,tol_golden);
        V1(a_c,z_c)     = fun_rhs(optimal_aprime);
        pol_ap(a_c,z_c) = optimal_aprime;
        
    end
end

end %end function "sub_vfi_onestep"

% function F = fun_rhs(aprime)
% 
% cons = cash-aprime;
% if cons<0
%     disp('cons<0')
%     keyboard
% end
% EV_interp = myinterp1(a_grid,EVz,aprime,1);
% F = util(cons,gamma) + beta*EV_interp;
% 
% end %end function "fun_rhs"

