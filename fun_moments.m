function [mom] = fun_moments(sol,agg,par)

z_grid = par.z_grid;
a_grid = par.a_grid;
mu     = sol.mu;
pol_c  = sol.pol_c;
ave_z  = par.Ls;

std_z = sqrt(sum((z_grid'-ave_z).^2.*mu,'all'));
cv_z  = std_z/ave_z;

ave_c = agg.CC;
std_c = sqrt(sum((pol_c-ave_c).^2.*mu,'all'));
cv_c  = std_c/ave_c;

ave_a = agg.Ks;
std_a = sqrt(sum((a_grid-ave_a).^2.*mu,'all'));
cv_a  = std_a/ave_a;

% Pack output
mom = v2struct(cv_z,cv_c,cv_a);

end

