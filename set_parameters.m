function par = set_parameters()

% Calibration Pontus Rendall

%% Set flags
do_GE     = 1; % Flag 0/1
do_save   = 1;
do_howard = 1;
disp_vfi  = 0;
disp_mu   = 0;
verbose   = 0;

%% Set economic parameters
alpha = 0.33;
beta  = 1.03^(-1/4);
delta = 0.025;  % Capital depreciation rate
gamma = 3;
repl  = 0.4; % Replacement rate for unemployment benefits
nss   = 0.95;
tau   = repl*(1-nss)/(nss+repl*(1-nss)); % Tax rate on labor income
d     = 0.1; % Job destruction rate
f     = d*nss/(1-nss+d*nss); % Job finding rate

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
a_max   = 400;
a_space = 3;
a_grid  = make_grid(a_min,a_max,na,a_space,1);

%% Exogenous shock z={e,u}
% - Grid:
    z_grid = [1,0]';
    nz = length(z_grid);
% - Transition matrix:
    pi_z  = [1-d d;f,1-f]; 
%   where   d is the separation rate, f is the job finding rate

n_ss = f/(d+f);

%% Set government taxes

%par.tau = par.repl*(1-n_ss)/(n_ss+par.repl*(1-n_ss));

%% Pack parameters into struct par
par = v2struct(alpha,beta,delta,gamma,repl,tau,na,a_grid,nz,z_grid,pi_z,...
    tol_v,maxit_v,do_howard,do_save,disp_vfi,tol_dist,maxiter_dist,disp_mu,...
    alg_mu,tol_golden,alg_lin,alg_howard,tol_ge,maxit_ge,verbose,do_GE,n_how,...
    alg_GE,damp,alg_vfi,alg_ergo);


end %end function "set_parameters"

