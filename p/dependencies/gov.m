classdef gov
    methods(Static)

        function net = tax(gross, lambda, tau)
            net = lambda*(gross.^(1-tau));
        end

        function Pr = getVoteShare(votes_today, adistr, amu, agrid)

            acond = compute.condense(adistr,amu,agrid);
            Pr = sum(sum(votes_today.*acond));

        end

        function R = getR(phis, Theta, r, prdraw)
            prob1 = 1 ./ (1 + exp(-phis(r) * (Theta - 0.5)));

            R = 2 - (prob1>=prdraw);

        end

    end
end