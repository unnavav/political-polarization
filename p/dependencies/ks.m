classdef ks
    methods(Static)

        function [Kprdata, Rdata, Prdata, V, G, EV] = getRegData(Rt, terms, vTol, verbose)

            T = length(Rt);
            pil = terms.pil;
            agrid = terms.agrid;
            lgrid = terms.lgrid;
            Kgrid = terms.Kgrid;
            phis = terms.phis;
            na = length(agrid); nl = length(lgrid); nm = length(Kgrid);

            Kprdata = zeros(T,1);
            Rdata = zeros(T,1);
            Prdata = zeros(T,1);
%             gamma = terms.gamma;

            nmu = na*10;
            amu = linspace(agrid(1), agrid(na), nmu);
    
            % solve for the household problem
            fprintf("Solving HH problem...\n")
            terms.lambda = .6; %from earlier versions of solving for SS
            [V, G, EV] = ks.solve(terms, vTol, verbose);

            % get inital conditions for the regression data
            Kprdata(1) = median(Kgrid);
            Votes = EV(:,1,:,:) >= EV(:,2,:,:);
            Votes = squeeze(Votes);
            g_int = squeeze(G(median(1:nm),1,:,:)); % TODO, figure out actual start
            [adistr, Kagg] = HH.getDist(g_int, amu, agrid, pil, false); 
            todays_voting_rule = squeeze(Votes(median(1:nm),:,:));
            Pr = gov.getVoteShare(todays_voting_rule, adistr, amu, agrid);
            Prdata(1) = Pr;
            Rdata(1) = gov.getR(phis,Pr,1,terms.rnseed);

            rng(terms.rnseed);

            fprintf("Generating regression data...\n")
            for t = 2:1:T
                K_t = Kprdata(t-1);
                pdraw = rand(1);
                R = Rdata(t-1); %voting outcome is decided by previous period
                % get the distribution

                [im, we] = compute.weight(Kgrid, K_t);

                g_int = we*G(im, R, :, :) + (1-we)*G(im+1, R, :, :);
                g_int = squeeze(g_int);
                todays_voting_rule = we*Votes(im,:,:) + ...
                    (1-we)*Votes(im+1, :, :);
                todays_voting_rule = squeeze(todays_voting_rule);


                [adistr, Kagg] = HH.getDist(g_int, amu, agrid, pil, false);

                Kprdata(t) = Kagg;
                Pr = gov.getVoteShare(todays_voting_rule, adistr, amu, agrid);
                Prdata(t) = Pr;
                Rdata(t) = gov.getR(phis,Pr,R,terms.rnseed);

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
            pil = terms.pil;
            g = terms.G;

            Kfore = terms.Kfore;
            Rfore = terms.Rfore;
            etagrid = terms.etagrid;
            taugrid = terms.taugrid;
            delta = terms.delta;
            captax = terms.captax;

            V = zeros(nm, nr, nl, na);
            G = V;
            V0 = V;
            TV = V; TG = G;

            % set up V so that it doesn't start empty
            % using believable r and w values: 4% interest, w = 1.3
            scale = .25;
            for ir = 1:nr 
                for im = 1:nm
                    for ia = 1:na
                        kval = agrid(ia);
                        for il = 1:nl
                            yval = scale*(1+0.04*(1-captax(il)))*kval + ...
                                1.3*lgrid(il) - 0.04*phi;
                            ymin = max(1e-10, yval);
                            V(im, ir, il, ia) = HH.u(ymin, sigma);
                        end
                    end
                end
            end

            dist = 1e5;
            iter_ct = 1;
            
            % forecasting K, R. outputs nkx1 and nkxr forecasts
            Kpr = ks.forecastK(Kfore, Kgrid); % nk x r
            Rpr = ks.forecastR(Rfore, Kgrid); % nk x r

            %then getting forecasted prices once, given future K.
            % nmxnr grids
            rgrid1 = vaas.calcr(alpha, delta, Kgrid, etagrid(1));
            rgrid2 = vaas.calcr(alpha, delta, Kgrid, etagrid(2));
            wgrid1 = vaas.calcw(alpha, Kgrid, etagrid(1));
            wgrid2 = vaas.calcw(alpha, Kgrid, etagrid(2));

            rgrid = [rgrid1' rgrid2']; wgrid = [wgrid1' wgrid2']; 

            while dist > vTol
                
                EV = ks.getExpectation(V, pil, Kpr, Rpr, Kgrid);

                % now converging on the value function and decision rule
                % for every capital-regime combo (EGM bc this is 50
                % convergences)

                for im = 1:nm
                    for ir = 1:nr
                        pol_terms = terms;
                        pol_terms.eta = terms.etagrid(ir);
                        pol_terms.tau = taugrid(ir);
                        pol_terms.r = rgrid(im, ir);
                        pol_terms.w = wgrid(im, ir);
                        [TV(im, ir, :,:), TG(im, ir,:,:)]= egm.solve(...
                            pol_terms, ...
                            squeeze(EV(im, ir, :,:)), ...
                            squeeze(V(im, ir,:,:)));
                    end 
                end

                % check distance
                dist = compute.dist(V, TV, 4);
                kdist = compute.dist(G, TG, 4);
            
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
        % moment of future capital. 
        function EV = getExpectation(V, pil, Kpr, Rpr, Kgrid)

            [nm, nr, ne, na] = size(V);
            EV1 = zeros(size(V));
            EV2 = EV1;
            EV3 = EV1;

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
            

            % step 2: weighting based on transition probability
            for im = 1:nm
                for ir = 1:nr
                    pr_R1 = Rpr(im, ir);
                    EV2(im, ir, :, :) = pr_R1*EV1(im, 1, :, :) + (1-pr_R1)*EV1(im, 2, :, :);
                end
            end

            % step 3: weighting based on forecasted K
            for im = 1:nm
                for ir = 1:nr
                    ix_m = ix(im, ir); we_m = we(im, ir); 
                    EV3(im, ir,:,:) = we_m*EV2(ix_m, ir, :, :) ...
                        + (1-we_m)*EV2(ix_m+1, ir, :, :);
                end
            end

            EV = EV3;
        end

        function kpr = forecastK(fore,k)            

            % preds = nk x 2
            preds = [ones(length(k),1), log(k)'];
            preds = preds'; % 2 x nk            

            preds1 = fore(1,:)*preds; % regime 1
            preds2 = fore(2,:)*preds; % regime 2
            preds = [preds1' preds2'];
            kpr = exp(preds);
        end

        function pr = forecastR(fore,k)

            % preds = nk x 2, fore = 2x2
            % there is a way to do this with matrix algebra that i cba to
            % figure out
            preds = [ones(length(k),1), log(k)'];
            pr1 = (1+exp(sum(fore(1,:).* preds, 2))).^-1;    
            pr2 = (1+exp(sum(fore(2,:).* preds, 2))).^-1;    
            pr = [pr1 pr2];

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