classdef gov
    methods(Static)

        function net = tax(gross, lambda, tau)
            net = gross - lambda*(gross.^(1-tau));
        end

    end
end