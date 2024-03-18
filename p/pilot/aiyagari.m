classdef aiyagari
    methods(Static)

        %% compute prices for given K and L
        function [r, w] = prices(K, L, alpha, delta)
           w = (1-alpha)*(K^alpha)*(L^(-alpha));
            r = (alpha)*((K/L)^(alpha-1)) - delta;
        end

        function u = uc(kpr, params)
            y = params(1);
            beta = params(2);
            sigma = params(3);
            ev = params(4);
            
            if sigma == 1
                u = log(y-kpr) + beta*ev;
            else
                u = (y- kpr) + beta*ev;
            end
        end

    end
end