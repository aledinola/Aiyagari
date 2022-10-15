function [] = make_plots_FK(sol,par)

mu_a = sum(sol.mu,2);
nz = par.nz;

a_lim = min(30,par.a_grid(end));
a_lim1 = min(30,par.a_grid(end));
a_lim2 = min(0.5,par.a_grid(end));
[~,loc] = min(abs(par.a_grid-a_lim));
[~,loc1] = min(abs(par.a_grid-a_lim1));
[~,loc2] = min(abs(par.a_grid-a_lim2));

figure
plot(par.a_grid,sol.pol_c(:,1),'linewidth',2)
hold on
plot(par.a_grid,sol.pol_c(:,nz),'linewidth',2)
legend('low z','high z')
xlabel('Assets')
ylabel('Consumption')
title('Policy function')
if par.do_save==1; print(fullfile('figs','pol_c1'),'-dpng'); end

figure
plot(par.a_grid(1:loc),sol.pol_c(1:loc,1),'linewidth',2)
hold on
plot(par.a_grid(1:loc),sol.pol_c(1:loc,2),'linewidth',2)
xlim([par.a_grid(1),par.a_grid(loc)])
legend('low z','high z','location','best')
xlabel('Assets')
ylabel('Consumption')
title('Policy function for consumption')
if par.do_save==1; print(fullfile('figs','pol_c'),'-dpng'); end

figure
plot(par.a_grid(1:loc),par.a_grid(1:loc),'--','linewidth',2)
hold on
plot(par.a_grid(1:loc),sol.pol_ap(1:loc,1),'linewidth',2)
hold on
plot(par.a_grid(1:loc),sol.pol_ap(1:loc,2),'linewidth',2)
xlim([par.a_grid(1),par.a_grid(loc)])
legend('45 line','low z','high z','location','best')
xlabel('Assets')
ylabel('Next-period assets')
title('Policy function for asset holdings')
if par.do_save==1; print(fullfile('figs','pol_ap'),'-dpng'); end

figure
plot(par.a_grid,sol.pol_ap(:,1),'linewidth',2)
hold on
plot(par.a_grid,sol.pol_ap(:,2),'linewidth',2)
hold on
plot(par.a_grid,sol.pol_ap(:,3),'linewidth',2)
hold on
plot(par.a_grid,sol.pol_ap(:,4),'linewidth',2)
xlim([par.a_grid(1),par.a_grid(loc2)])
legend('z1','z2','z3','z4','location','northwest')
xlabel('Wealth a')
ylabel('Savings for tomorrow a'' ')
title('Policy function for asset holdings')
if par.do_save==1; print(fullfile('figs','pol_ap2'),'-dpng'); end



figure
plot(par.a_grid,sol.mu(:,1),'linewidth',2)
hold on
plot(par.a_grid,sol.mu(:,2),'linewidth',2)
xlim([par.a_grid(1),par.a_grid(loc1)])
legend('low z','high z','location','best')
xlabel('Assets')
ylabel('Density')
title('Distribution')
if par.do_save==1; print(fullfile('figs','mu'),'-dpng'); end

figure
plot(par.a_grid,mu_a,'linewidth',2)
xlabel('Assets')
ylabel('Density')
title('Distribution')
if par.do_save==1; print(fullfile('figs','mu_a1'),'-dpng'); end

figure
plot(par.a_grid,mu_a,'linewidth',2)
xlim([par.a_grid(1),par.a_grid(loc1)])
xlabel('Assets')
ylabel('Density')
title('Distribution')
if par.do_save==1; print(fullfile('figs','mu_a'),'-dpng'); end

end %end function

