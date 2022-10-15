function [V1,pol_ap] = sub_vfi_onestepx(V0,cash_arr,par)

% DESCRIPTION
%   sub_vfi_onestep implements the Bellman operator: given an initial V0 it
%   performs the maximization and returns an updated V1. We use goldenx
%   instead of golden. Compare to sub_vfi_onestep.
% INPUTS
%   V0         : Value function (initial condition), dim: (na,nz)
%   cash_arr   : Precomputed cash-on-hand, dim: (na,nz)
%   par        : Struct with model parameters
% OUTPUTS
%   V1         : Updated value function, dim: (na,nz)
% AUXILIARY
%   goldenx    : this function does vectorized golden search method.
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

for z_c = 1:nz
    EVz = EV(:,z_c);
    %for a_c = 1:na
        
        % Set lower and upper bound for a'. The upper bound cannot be
        % larger than max(a_grid); so we avoid doing extrapolation (risky!)
        cash_vec  = cash_arr(:,z_c);
        ap_lb_vec = repmat(ap_lb,na,1);
        ap_ub_vec = min(ap_ub,cash_vec-1e-10);
        
        fun_rhs = @(ap) fun.util(cash_vec-ap,gamma) + beta*interp1qr(a_grid,EVz,ap);
        
        % Maximize the RHS of the Bellman equation
        optimal_aprime  = goldenx(fun_rhs,ap_lb_vec,ap_ub_vec,tol_golden);
        V1(:,z_c)     = fun_rhs(optimal_aprime);
        pol_ap(:,z_c) = optimal_aprime;
        
    %end
end

end %end function "sub_vfi_onestepx"

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

