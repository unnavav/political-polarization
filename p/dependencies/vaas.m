classdef vaas
    methods(Static)

        function r = calcr(alpha, delta, k, eta)
            r = (alpha)*(k/(1+eta))^(alpha - 1) - delta;
        end

        function w = calcw(alpha, k, eta)
            w = (1-alpha)*(k/(1+eta))^alpha;
        end

    end
end