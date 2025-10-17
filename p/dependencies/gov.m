classdef gov
    methods(Static)

        function net = tax(gross, lambda, tau)
            net = lambda*(gross.^(1-tau));
        end

        function Pr = getVoteShare(Vote_KR, adistr, amu, agrid)

            acond = compute.condense(adistr,amu,agrid);
            Pr = sum(sum(Vote_KR.*acond));

        end

        function R = getR(phis, Theta, r, rnseed)
            
            rng(rnseed);
            prob1 = (exp(-phis(r)*(Theta-0.5)))^-1;
            draw = rand(1);

            R = 2 - (prob1>=draw);

        end

    end
end