%% Aiayagari model
% Please set model parameters (economic parameters, flag, computational
% parameters, grids, etc.) in the function set_parameters.m
% All model parameters are collected in the struct "par".
% I provide two different calibrations for the Aiyagari model:
% - set_parameters_FK: used in chapter 9 of Fehr and Kindermann
% - set_parameters_ACI: used in Acikgoz JET (2018) to show multiplicity of
%   equilibria in a Bewley model with production. This calibration uses
%   quite extreme parameter value, e.g. capital depreciation is 100%.

clear
clc
close all
addpath(genpath(fullfile('tools')));
format long g

disp("Hello world!")

%% Set parameters

%par = set_parameters();
par = set_parameters_FK();
%par = set_parameters_ACI();

% Initial guess for the interest rate
rl = 0;
rh = min((1/par.beta-1)*0.98,8);
r0 = 0.040124246749556;

err_ge  = 1+par.tol_ge;
iter_ge = 0;

options = optimset('TolX',par.tol_ge);

if par.do_GE==0
    disp('PARTIAL EQUILIBRIUM')
    %tic
    [ED,sol,agg,par] = excess_demand(r0,par);
    %toc
elseif par.do_GE==1
    disp('GENERAL EQUILIBRIUM')
    tic
    switch par.alg_GE
        case 1
            disp('...with Brent method (fzero)')
            r = fzero(@(x) excess_demand(x,par),[rl,rh],options);
            %r = fzero(@(x) excess_demand(x,par),r0,options);
        case 2
            disp('...with Illinois method')
            r = illinois(@(x) excess_demand(x,par),rl,rh,par.tol_ge,1e-6,par.maxit_ge);
        case 3
            % Dampening or Gauss Seidel
            iter_ge = 0;
            err_ge  = 1+par.tol_ge;
            r = r0;
            while err_ge>par.tol_ge && iter_ge<=par.maxit_ge
                
                r_implied = excess_demand(r,par);
                r = par.damp*r+(1-par.damp)*r_implied;
                err_ge = abs(r_implied-r);
                iter_ge = iter_ge+1;
                
            end
            
        case 4
            % GRAPHICAL ANALYSIS
            disp('GRAPHICAL ANALYSIS')
            % Draw capital demand and capital supply as functions of r
            nr = 30;
            r_vec  = linspace(0,8,nr)';
            Kd_vec = zeros(nr,1);
            Ks_vec = zeros(nr,1);
            
            for ir = 1:nr
                fprintf('Solve the model at r = %f \n',r_vec(ir))
                [~,~,~,~,Kd,Ks] =  excess_demand(r_vec(ir),par);
                Kd_vec(ir) = Kd;
                Ks_vec(ir) = Ks;
            end
            
            [~,r_loc] = min(abs(Kd_vec-Ks_vec));
            r = r_vec(r_loc);
            
            figure
            plot(Kd_vec,r_vec,'k','linewidth',2)
            hold on
            plot(Ks_vec,r_vec,'--k','linewidth',2)
            legend('K(r), capital demand','A(r), capital supply')
            axis tight, grid on
            xlabel('Capital','fontsize',14)
            ylabel('Interest rate r','fontsize',14)
            hold off
            if par.do_save==1; print(fullfile('figs','multi_eq'),'-dpng'); end
            
        otherwise
            error('alg_GE is not specified correctly!')
    end
    toc
    disp('******* GE found! *****************')
    [ED,sol,agg,par] = excess_demand(r,par);
end
% Compute second moments: coef variation, Gini, etc.
[mom] = fun_moments(sol,agg,par);

write_output(ED,agg,mom,par);

make_plots_FK(sol,par);
