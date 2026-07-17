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

        function [EV_R1, EV_R2, Votes_EV] = getVotingExpectations(V, pil, pid, Kpr, Kgrid)
            [nm, nr, nd, ne, na] = size(V);
            EV = zeros(size(V));

            [ix, we] = ks.weight(Kgrid,Kpr);

            % step 1: updating EV(a,e)
            for id = 1:nd
                for im = 1:nm
                    for ir = 1:nr
                        for ia = 1:na
                            for ie = 1:ne
                                EV(im, ir, id, ie, ia) = pil(ie,:)*squeeze(V(im, ir, id, :, ia));
                            end
                        end
                    end
                end
            end

            % step 2: updating with pid
            for id = 1:nd
                EV(:,:,id,:,:) = pid(id,1)*EV(:,:,1,:,:) + ...
                    (1-pid(id, 1))*EV(:,:,2,:,:);
            end
            

            EV_R1 = zeros(size(V));  % value if NEXT regime is forced to 1
            EV_R2 = zeros(size(V));  % value if NEXT regime is forced to 2
        
            for id = 1:nd
                for im = 1:nm
                    for ir = 1:nr
                        ix_m = ix(id, im, ir); 
                        we_m = we(id, im, ir);
            
                        % Interpolate from the ε-averaged EV *holding R' fixed*
                        EV_R1(im, ir, id, :, :) = ...
                            we_m    * EV(ix_m, 1, id, :, :) + ...
                            (1-we_m)* EV(ix_m+1, 1, id, :, :);
            
                        EV_R2(im, ir, id, :, :) = ...
                            we_m    * EV(ix_m,   2, id, :, :) + ...
                            (1-we_m)* EV(ix_m+1, 2, id, :, :);
                    end
                end
            end
        
            % Step 3: voting rule (elementwise comparison over (ε,a))
            Votes_EV = (EV_R1 >= EV_R2);    % nm x nr x nl x na

        end


    end
end