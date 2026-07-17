classdef vaas
    methods(Static)

        function r = calcr(alpha, delta, k, eta)
            r = (alpha).*(k./(1+eta)).^(alpha - 1) - delta;
        end

        function w = calcw(alpha, k, eta)
            w = (1-alpha).*(k./(1+eta)).^alpha;
        end

        function dw = dwde(alpha, k, eta)
            dw = -alpha*(1-alpha).*(k).^alpha.*(1+eta).^(-alpha - 1);
        end

        function dr = drde(alpha, k, eta)
            dr = alpha*(1-alpha).*(k).^(alpha-1).*(1+eta).^(-alpha - 1);
        end
    end
end