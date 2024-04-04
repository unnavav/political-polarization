classdef aiyagari
    methods(Static)

        %% compute prices for given K and L
        function [r, w] = prices(K, L, alpha, delta)
           w = (1-alpha)*(K^alpha)*(L^(-alpha));
            r = (alpha)*((K/L)^(alpha-1)) - delta;
        end

        function u = uc(apr, params)
            y = params(1);
            beta = params(2);
            sigma = params(3);
            
            cons = y-apr;
            if cons <= 0
                u = 2000;
            elseif sigma == 1
                u = log(cons);
            else
                u = (cons)^(1-sigma)/(1-sigma);
            end
        end

    end
end