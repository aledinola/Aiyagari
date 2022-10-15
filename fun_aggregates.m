function [agg] = fun_aggregates(sol,par)

% UNPACK
mu     = sol.mu;     % Distribution mu, dim: (na,nz)
a_grid = par.a_grid; % Assets grid,     dim: (na,1)
%z_grid = par.z_grid; % Labor prod grid  dim: (nz,1)
pol_c  = sol.pol_c;  % Policy consump,  dim: (na,nz)
alpha  = par.alpha;
delta  = par.delta;

% Capital supply is equal to total savings
Ks = sum(a_grid.*mu,'all');

% Labor supply is exogenous. Check this in Pontus model: it is exogenous
% but not necessarily equal to one!
Ls = par.Ls;

% Aggregate output
YY = Ks^alpha*Ls^(1-alpha);

% Aggregate consumption
CC = sum(pol_c.*mu,'all');

% Aggregate investment
II = delta*Ks;

% Market clearing residual
LHS = YY;
RHS = CC+II;
res = abs(LHS-RHS);
res1 = abs(YY-par.wage*Ls-(par.r+delta)*Ks);

% Pack outputs
agg.Ks = Ks;
agg.Ls = Ls;
agg.YY = YY;
agg.CC = CC;
agg.II = II;
agg.res = res;
agg.res1 = res1;

end %end function fun_aggregates

