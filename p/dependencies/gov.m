classdef gov
    methods(Static)

        function net = tax(gross, lambda, tau)
            net = lambda*(gross.^(1-tau));
        end

    end
end