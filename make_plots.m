function [] = make_plots(sol,par)

a_lim = min(30,par.a_grid(end));
a_lim1 = min(130,par.a_grid(end));
[~,loc] = min(abs(par.a_grid-a_lim));
[~,loc1] = min(abs(par.a_grid-a_lim1));

figure
plot(par.a_grid,sol.pol_c(:,1),'linewidth',2)
hold on
plot(par.a_grid,sol.pol_c(:,2),'linewidth',2)
legend('Employed','Unemployed')
xlabel('Assets')
ylabel('Consumption')
title('Policy function')
if par.do_save==1; print(fullfile('figs','pol_c1'),'-dpng'); end

figure
plot(par.a_grid(1:loc),sol.pol_c(1:loc,1),'linewidth',2)
hold on
plot(par.a_grid(1:loc),sol.pol_c(1:loc,2),'linewidth',2)
xlim([par.a_grid(1),par.a_grid(loc)])
legend('Employed','Unemployed','location','best')
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
legend('45 line','Employed','Unemployed','location','best')
xlabel('Assets')
ylabel('Next-period assets')
title('Policy function for asset holdings')
if par.do_save==1; print(fullfile('figs','pol_ap'),'-dpng'); end

figure
plot(par.a_grid(1:loc1),sol.mu(1:loc1,1),'linewidth',2)
hold on
plot(par.a_grid(1:loc1),sol.mu(1:loc1,2),'linewidth',2)
xlim([par.a_grid(1),par.a_grid(loc1)])
legend('Employed','Unemployed','location','best')
xlabel('Assets')
ylabel('Density')
title('Distribution')
if par.do_save==1; print(fullfile('figs','mu'),'-dpng'); end

end %end function

