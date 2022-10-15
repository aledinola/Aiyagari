function [F,cons]  = ReturnFn(a_val,z_val,aprime_val,par)

gamma = par.gamma;
cash = fun.fun_cash(a_val,z_val,par); % dim: (a,z)
cons = cash-aprime_val; % dim: (a,z,a')

F = repmat(-inf,size(cons));
idc = cons>0;
F(idc) = fun.util(cons(idc),gamma);

end