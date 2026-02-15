classdef ks
    methods(Static)

        function [Kprdata, Rdata, Prdata, ddata, distr_array, V, G, EV] = ...
                getRegData(T, terms, vTol, verbose)

            agrid = terms.agrid;
            lgrid = terms.lgrid;
            Kgrid = terms.Kgrid;
            na = length(agrid); nl = length(lgrid); nm = length(Kgrid);

            pil = terms.pil;
            pid = terms.pid;

            % forecasting K, R. outputs nkx1 and nkxr forecasts
            Kpr = ks.forecastK(terms.Kfore, Kgrid); % nk x r
            Rpr = ks.forecastR(terms.Rfore, Kgrid); % nk x r
            terms.Kpr = Kpr;
            terms.Rpr = Rpr;

            nmu = na*10;
            amu = linspace(agrid(1), agrid(na), nmu);
    
            % solve for the household problem
            fprintf("Solving HH problem...\n")
            [V, G, EV] = ks.solve(terms, vTol, verbose);

            % Votes  
            [EV1, EV2, Votes_EV] = gov.getVotingExpectations(V, pil, pid, Kpr, Kgrid);

            % get inital conditions for the regression data
            Kprdata = zeros(T,1);
            distr_array = cell(T,1);
            g0 = terms.starter_distr;
            g0_cond = compute.condense(g0, amu, agrid);
            K0 = dot(squeeze(sum(sum(g0_cond,1),2)),agrid);
            Kprdata(1) = K0; distr_array{1} = g0;

            % start with regime one, then check
            Rdata = ones(T,1);
            Prdata = Rdata;

            ddata = predict.sim(T,2,"default",terms.pid);

            fprintf("Generating regression data...\n")
            for t = 2:1:T
                Kt = Kprdata(t-1);
                Rt = Rdata(t);
                dt = ddata(t);
                g_prev = distr_array{t-1};
                
                [ix, we] = compute.weight(Kgrid, Kt);
                g_t = we*G(ix, Rt, dt, :,:) + (1-we)*G(ix+1, Rt, dt, :, :);
                g_t = squeeze(g_t);
                g_today = HH.transitDistr(g_t, g_prev, amu, agrid, pil);

                distr_array{t} = g_today;
                acond = compute.condense(g_today, amu, agrid);
                Kpr = dot(squeeze(sum(sum(acond,1),2)),agrid);
                Kprdata(t) = Kpr;

                % now use K today, Kpr, and R today to back out max vote,
                % along with the actual distribution over wealth. Note that
                % this voting rule already interpolates over the future
                % forecast, so I only need to interpolate over today's K,
                % then force it back to binary (otherwise there's a decimal
                % value on whether or not I'll vote for R = 1)

                todays_votes = we*Votes_EV(ix, Rt, dt, :, :) + ...
                    (1-we)*Votes_EV(ix+1, Rt, dt, :, :);
                todays_votes = squeeze(todays_votes);
                todays_votes = (todays_votes >= .5);

                vote_total = sum(sum(sum(squeeze(acond(dt,:,:)).*todays_votes)))*2;
                %multiplying by 2 bc each dimension has 50% of mass
                Prdata(t) = vote_total;
                if (vote_total <=.5) 
                    Rdata(t+1) = 2; 
                else 
                    Rdata(t+1) = 1;
                end

                if mod(t,100) == 0
                    fprintf("\n\t t = %i", t)
                end
            end
        end


        function [V, G, V0] = solve(terms, vTol, verbose)

            nr = length(terms.etagrid);
            alpha = terms.alpha;
            sigma = terms.sigma;
            phi = terms.phi;

            lgrid = terms.lgrid; nl = length(lgrid);
            agrid = terms.agrid; na = length(agrid);
            Kgrid = terms.Kgrid; nm = length(Kgrid);
            dgrid = terms.dgrid; nd = length(dgrid);
            pil = terms.pil; pid = terms.pid;
            g = terms.G;

            etagrid = terms.etagrid;
            taugrid = terms.taugrid;
            captax = terms.captax;

            V = zeros(nm, nr, nd, nl, na);
            G = V;
            V0 = V;
            TV = V; TG = G;

            % set up V so that it doesn't start empty
            % using believable r and w values: 4% interest, w = 1.3
            scale = .25;
            for id = 1:nd
                for ir = 1:nr 
                    for im = 1:nm
                        for ia = 1:na
                            kval = agrid(ia);
                            for il = 1:nl
                                yval = scale*(1+0.04*(1-captax(il)))*kval + ...
                                    1.3*lgrid(il) - 0.04*phi;
                                ymin = max(1e-10, yval);
                                V(im, ir, id, il, ia) = HH.u(ymin, sigma);
                            end
                        end
                    end
                end
            end

            dist = 1e5;
            iter_ct = 1;
            

            %then getting forecasted prices once, given future K.
            % nmxnr grids
            rgrid = zeros(nm, nr, nd);
            wgrid = rgrid;
            for ir = 1:nr
                for id = 1:nd
                   rgrid(:,ir,id)  = vaas.calcr(alpha, dgrid(id), Kgrid, etagrid(ir));
                   wgrid(:,ir,id) = vaas.calcw(alpha, Kgrid, etagrid(ir));
                end
            end

            % calculating lambda from prices.
            lambda_grid  = (wgrid(:,:,1) .^ terms.taugrid) .* ...
                (ones(size(wgrid,1),1) * terms.lamval);

            Kpr = terms.Kpr;
            Rpr = terms.Rpr;

            while dist > vTol
                
                EV = ks.getExpectation(V, pil, pid, Kpr, Rpr, Kgrid);

                % now converging on the value function and decision rule
                % for every capital-regime combo (EGM bc this is 50
                % convergences)

                for id = 1:nd
                    for im = 1:nm
                        for ir = 1:nr
                            pol_terms = terms;
                            pol_terms.eta = terms.etagrid(ir);
                            pol_terms.tau = taugrid(ir);
                            pol_terms.r = rgrid(im, ir, id);
                            pol_terms.w = wgrid(im, ir, id);
                            pol_terms.lamval = lambda_grid(im, ir);
                            [TV(im, ir, id, :,:), TG(im, ir, id,:,:)]= ...
                                egm.solve(pol_terms, ...
                                squeeze(EV(im, ir, id, :,:)), ...
                                squeeze(V(im, ir, id,:,:)));
                        end 
                    end
                end

                % check distance
                dist = compute.dist(V, TV, 5);
                kdist = compute.dist(G, TG, 5);
            
                if mod(iter_ct, 25) == 0
                    fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
                        "\n\t\t||TG - G|| = %4.6f", iter_ct, dist, kdist);
% %                     fprintf("\nInitial Values:");
% %                     fprintf("\nMin Kgrid: %2.4f, Max Kgrid: %2.4f", min(Kgrid), max(Kgrid));
% %                     fprintf("\nInitial Capital Forecasts: Kp = %2.4f, Kl = %2.4f", Kp, Kl);
% %                     fprintf("\nInitial Policy Function: Min Gp = %2.4f, Max Gp = %2.4f", ...
% %                         min(TGP(:)), max(TGP(:)));
% %                      fprintf("\nInitial Policy Function: Min Gl = %2.4f, Max Gl = %2.4f", ...
% %                         min(TGL(:)), max(TGL(:)));
% %                    fprintf("\nInitial Value Function: Min VP = %2.4f, Max VP = %2.4f", ...
% %                        min(TVP(:)), max(TVP(:)));
% %                    fprintf("\nInitial Value Function: Min VL = %2.4f, Max VL = %2.4f", ...
% %                        min(TVL(:)), max(TVL(:)));
                end
            
                iter_ct = iter_ct + 1;
            
                G = 0.4 * TG + 0.6 * G;
                V = 0.4 * TV + 0.6 * V;

            end

            if verbose
                fprintf("\n\tIteration %i: ||TV - V|| = %4.6f" + ...
                    "\t||TG - G|| = %4.6f\n", iter_ct, dist, kdist);
            end

        end

        % get expectation
        % V is a four dimensional object-- K x R x a x e--- and so when
        % we're forming expectations, there's three steps: predict v(a,e)
        % given the transition grid. Then weight between the two possible
        % regimes based on the probability of transition given current
        % regime & capital. Then finally interpolate across the expected
        % moments of future capital. And I added weighting on a capital
        % depreciation shock. Is it time for me to give up I want to give
        % up
        function EV = getExpectation(V, pil, pid, Kpr, Rpr, Kgrid)

            [nm, nr, nd, ne, na] = size(V);
            EV1 = zeros(size(V));
            EV2 = EV1;
            EV3 = EV1;
            EV4 = EV3;

            [ix, we] = ks.weight(Kgrid,Kpr);

            % step 1: updating EV(a,e)
            for id = 1:nd
                for im = 1:nm
                    for ir = 1:nr
                        for ia = 1:na
                            for ie = 1:ne
                                EV1(im, ir, id, ie, ia) = pil(ie,:)*squeeze(V(im, ir, id, :, ia));
                            end
                        end
                    end
                end
            end
            
            % step 2: weighting based on transition probability
            for id = 1:nd
                for im = 1:nm
                    for ir = 1:nr
                        p_next1 = Rpr(id,im,ir);            % P(R' = 1 | K,R,δ)
                        EV2(im, ir, id, :, :) = p_next1*EV1(im, 1, id, :, :) + ...
                            (1-p_next1)*EV1(im, 2, id, :, :);
                    end
                end
            end

            % step 3: weighting based on forecasted K
            for id = 1:nd
                for im = 1:nm
                    for ir = 1:nr
                        ix_m = ix(id, im, ir); we_m = we(id, im, ir); 
                        EV3(im, ir, id, :,:) = we_m*EV2(ix_m, ir, id, :, :) ...
                            + (1-we_m)*EV2(ix_m+1, ir, id, :, :);
                    end
                end
            end
            
            %step 4: weighting based on delta transition probability
            for id = 1:nd
                EV4(:,:,id,:,:) = pid(id,1)*EV3(:,:,1,:,:) + ...
                    (1-pid(id,1))*EV3(:,:,2,:,:);
            end

            EV = EV4;
        end

        function kpr = forecastK(fore,k)            

            [nd, ~, ~] = size(fore);
            kpr = zeros(nd, length(k), 2);

            for id = 1:nd
                % preds = nk x 2
                preds = [ones(length(k),1), log(k)'];
                preds = preds'; % 2 x nk            
    
                dfore = squeeze(fore(id, :, :));
                preds1 = dfore(1,:)*preds; % regime 1
                preds2 = dfore(2,:)*preds; % regime 2
                preds = [preds1' preds2'];
                kpr(id, :, :) = exp(preds);
            end
        end



        function pr = forecastR(fore,k)

            [nd, ~, ~] = size(fore);
            pr = zeros(nd, length(k), 2);

            % preds = nk x 2, fore = 2x2
            % there is a way to do this with matrix algebra that i cba to
            % figure out
            for id = 1:nd
                dfore = squeeze(fore(id, :, :));
                preds = [ones(length(k),1), log(k)'];
                z1 = preds * dfore(1,:)';            % nk x 1
                z2 = preds * dfore(2,:)';            % nk x 1
            
                pr1 = 1./(1 + exp(-z1));            % P(R' = 1 | R = 1, K)
                pr2 = 1./(1 + exp(-z2));            % P(R' = 1 | R = 2, K)
            
                pr(id, :, :)  = [pr1 pr2];                    % nk x 2, each entry in (0,1)
            end
        end

        function [ixmat, wemat] = weight(Kgrid, Kpr)

            Kpr = min(Kpr, ones(size(Kpr))*max(Kgrid));
            Kpr = max(Kpr, ones(size(Kpr))*min(Kgrid));
            nm = length(Kgrid);

            % Kpr = (nR*nd) x nm
            % creating index matrix first


            ixmat = discretize(Kpr, Kgrid);
            ixmat = min(ixmat, nm-1); % making sure it doesn't exceed

            Klower = Kgrid(ixmat); 
            Kupper = Kgrid(ixmat+1);

            wemat = (Kpr - Klower)./(Kupper - Klower);

        end

    end
end