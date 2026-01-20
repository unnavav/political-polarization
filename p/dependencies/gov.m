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

        function [EV1, EV2, Votes_EV] = getVotingExpectations(V, pil, Kpr, Kgrid)
            [nm, nr, ne, na] = size(V);
            EV1 = zeros(size(V));
            EV2 = EV1;

            [ix, we] = ks.weight(Kgrid,Kpr);

            % step 1: updating EV(a,e)
            for im = 1:nm
                for ir = 1:nr
                    for ia = 1:na
                        for ie = 1:ne
                            EV1(im, ir, ie, ia) = pil(ie,:)*squeeze(V(im, ir, :, ia));
                        end
                    end
                end
            end
            

            EV_R1 = zeros(size(V));  % value if NEXT regime is forced to 1
            EV_R2 = zeros(size(V));  % value if NEXT regime is forced to 2
        
            for im = 1:nm
                for ir = 1:nr
                    ix_m = ix(im, ir); 
                    we_m = we(im, ir);
        
                    % Interpolate from the ε-averaged EV *holding R' fixed*
                    EV_R1(im, ir, :, :) = ...
                        we_m    * EV1(ix_m,   1, :, :) + ...
                        (1-we_m)* EV1(ix_m+1, 1, :, :);
        
                    EV_R2(im, ir, :, :) = ...
                        we_m    * EV1(ix_m,   2, :, :) + ...
                        (1-we_m)* EV1(ix_m+1, 2, :, :);
                end
            end
        
            % Step 3: voting rule (elementwise comparison over (ε,a))
            Votes_EV = (EV_R1 >= EV_R2);    % nm x nr x nl x na

        end


    end
end