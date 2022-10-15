classdef fun
    %This class contains all functions related to the model
    % Written by Alessandro Di Nola
    % aledinola13@gmail.com
            
    methods (Static)
        
        function [U] = util(cons,gamma)
            % gamma is CRRA coefficient
            U = (cons.^(1-gamma)-1)/(1-gamma);
        end
        
        function [wage,KK] = prices(r,par)
            % Gives wage and aggregate capital from the interest rate
            alpha = par.alpha;
            delta  = par.delta;
            Ls     = par.Ls;
            
            KK = (alpha/(r+delta))^(1/(1-alpha))*Ls;
            wage = (1-alpha)*(KK/Ls)^alpha;
        end
        
        function [r] = intrate(KK,par)
            % Gives interest rate as a function of aggregate capital
            alpha = par.alpha;
            delta  = par.delta;
            Ls     = par.Ls;
            r = alpha*(KK/Ls)^(alpha-1)-delta;
        end
        
        function [cash] = fun_cash(a_val,z_val,par)
            % DESCRIPTION
            % fun_cash computes cash-on-hand defined from the household's budget
            % constraint as follows:
            %          c + a' = cash
            % INPUTS
            %   a_val      : Current assets, a, scalar.
            %   z_val      : Current exogenous shock, z, scalar.
            %   par        : Struct with model parameters.
            % OUTPUTS
            %   cash       : Precomputed cash-on-hand, scalar
            
            % Unpack
            r    = par.r;
            wage = par.wage;
            
            inc   = z_val*wage;
            cash  = (1+r)*a_val+inc;
            
        end %end function "fun_cash"
        
    end %end methods
end %end class

