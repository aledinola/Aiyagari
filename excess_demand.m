function [out,sol,agg,par,K_demand,K_supply] = excess_demand(r_input,par)

% DESCRIPTION
%   excess_demand computes the excess demand (or supply) for capital for a
%   given value of the interest rate.
% INPUTS
%   r_input :  Interest rate, scalar
%   par     :  Struct with parameters
% OUTPUTS
%   out     :  Excess demand, scalar
%   sol     :  Struct with value,policy functions and distribution
%   agg     :  Struct with scalar fields with model aggregates
%   par     :  Struct with parameters (fields for r and wage are updated)
%
% Written by Alessandro Di Nola
% aledinola13@gmail.com

verbose = par.verbose;

% Compute prices
par.r    = r_input;
[par.wage,K_demand] = fun.prices(r_input,par);

% Compute VFI
if verbose>=1
    disp('Start VFI...')
    tic
end
switch par.alg_vfi
    case 1
        disp('VFI: continuous optimization with interpolation')
        [sol,flag_vfi,iter_final] = compute_vfi(par);
    case 2 
        disp('VFI: discrete grid')
        [sol,flag_vfi,iter_final] = compute_vfi_discrete(par);
    otherwise
        error('alg_vfi is not specified correctly')
end
if verbose>=1
    disp('VFI running time:')
    toc
    fprintf('VFI terminated after %d iterations \n',iter_final)
end
if flag_vfi<0
    error('VFI did not converge!')
end

% Compute distribution
if verbose>=1
    disp('Start distribution...')
    tic
end
switch par.alg_mu
    case 1
        disp('Distribution algorithm: with loops')
        [mu,flag_mu,iter_final] = compute_mu(par.a_grid,sol.pol_ap,par.pi_z,par);
    case 2
        disp('Distribution algorithm: vectorized')
        %[mu,flag_mu,iter_final] = compute_mu_vec(par.a_grid,sol.pol_ap,par.pi_z,par);
        [mu,flag_mu,iter_final] = compute_mu_vec1(par.a_grid,sol.pol_ap,par.pi_z,par);
end

sol.mu = mu;
if verbose>=1
    toc
    fprintf('Distribution terminated after %d iterations \n',iter_final)
end
if flag_mu<0
    warning('Distribution did not converge!')
end

% Compute aggregates
agg = fun_aggregates(sol,par);

K_supply = agg.Ks; % capital supplied by households

if par.alg_GE==3
    r_implied = max(-0.01,min(fun.intrate(K_supply,par),0.8));
    out = r_implied;
    fprintf('-----------------------------------------------------------\n')
    %fprintf('r    = %f, r_implied = %f, err = %f \n',r,r_implied,out)
    fprintf('r    = %f, r_implied = %f, err = %f \n',r_input,r_implied,abs(r_input-r_implied))
    %fprintf('iter = %d, err_ge    = %f \n',iter_ge,err_ge)
    fprintf('-----------------------------------------------------------\n')
else
    out = (K_supply-K_demand)/K_supply;
    fprintf('-----------------------------------------------------------\n')
    %fprintf('r    = %f, r_implied = %f, err = %f \n',r,r_implied,out)
    fprintf('r    = %f, err = %f \n',r_input,out)
    %fprintf('iter = %d, err_ge    = %f \n',iter_ge,err_ge)
    fprintf('-----------------------------------------------------------\n')
end



end %end function "excess demand"

