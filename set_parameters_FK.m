function par = set_parameters_FK()

% Calibration Kindermann chapter 9

%% Set flags
do_GE     = 1; % Flag 0/1
do_save   = 1;
do_howard = 1;
disp_vfi  = 0;
disp_mu   = 0;
verbose   = 0;

%% Set economic parameters
gamma    = 2.0; % CRRA utility consumption
alpha    = 0.36;% capital coeff Cobb-Douglas
beta     = 0.96; % discount factor
delta    = 0.08;
rho_z    = 0.6;
sigma_z  = sqrt(0.04*(1-rho_z^2)); %st. deviation of AR1 shocks

%% Set computational parameters

% General equilibrium loop
tol_ge   = 1e-10;
maxit_ge = 300;
alg_GE   = 1;  % 1=Brent, 2=Illinois, 3=dampening, 4=draw capital demand and cap supply
damp     = 0.99;  % Weight on old iterate (used only if alg_GE=3)

% Value function
tol_v      = 1e-10;
maxit_v    = 3000;
tol_golden = 1e-6; 
alg_vfi    = 1; % 1=continuous optim with golden, 2=pure discrete grid
alg_lin    = 1; % Algorithm to solve large linear system in sub_howard_vec or in VFI discrete
alg_howard = 1; %1=with loops (sub_howard), 2=vectorized (sub_howard_vec)
n_how      = 100; %No. of Howards iterations (only if sub_howard is used)
% NOTE: with large state space, alg_howard=1 may be preferable.


% Distribution
tol_dist     = 1e-10;
maxiter_dist = 100000;
alg_mu       = 2; % Algorithm to solve for distrib: 1=with loops, 2=vectorized
alg_ergo     = 1; % Algorithm to find the ergodic distrib of Markov chain (only if alg_mu=2)

%% Grid for assets
na      = 1000; % No. of points for assets
a_min   = 0;
a_max   = 100;
a_space = 3;%0.01;
a_method = 1; % 1=nonlinspace, 2=Kindermann
a_grid  = make_grid(a_min,a_max,na,a_space,a_method);

%% Exogenous shock z
nz = 19;
ave_z = 0;
[z_grid_logs,pi_z] = rouwenhorst(nz,ave_z,rho_z,sigma_z);
z_grid = exp(z_grid_logs);
temp = pi_z^1000;
z_prob = temp(1,:)'; % (nz,1) vector with ergodic distrib of z
Ls = sum(z_grid.*z_prob); % average labor supply, fixed since it is inelastic

%% Set government taxes

%par.tau = par.repl*(1-n_ss)/(n_ss+par.repl*(1-n_ss));

%% Pack parameters into struct par
par = v2struct(gamma,alpha,beta,delta,Ls,na,a_grid,nz,z_grid,pi_z,...
    tol_v,maxit_v,do_howard,do_save,disp_vfi,tol_dist,maxiter_dist,disp_mu,...
    alg_mu,tol_golden,alg_lin,alg_howard,tol_ge,maxit_ge,verbose,do_GE,n_how,...
    alg_ergo,alg_GE,damp,alg_vfi);


end %end function "set_parameters"

