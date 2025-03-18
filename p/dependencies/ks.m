classdef ks
    methods(Static)

        function [Kdata, adist] = getRegData(jt, terms, foreguess, nm, dTol, vTol, verbose)

            T = length(jt);
            agrid = terms.agrid;
            lgrid = terms.lgrid;
            na = length(agrid); nl = length(lgrid);

            Kdata = zeros(1,T);

            gamma = terms.gamma;

            nmu = na*10;
            amu = linspace(agrid(1), agrid(na), nmu);

            terms.lamval = 0.2;
            [~, G, ~, EV] = ks.solve(nm, nl, na, terms, vTol, verbose);

            for t = 1:1:T

                K = jt(t);
                % get the distribution

                [adistr, Kagg] = HH.getDist(G, amu, agrid, pil, false);
                Kdata = Kagg;

            end

        end

        function [V, G, C, V0] = solve(nm, nl, na, terms, vTol, verbose)

            beta = terms.beta;
            sigma = terms.sigma;
            phi = terms.phi;
            lgrid = terms.lgrid;
            agrid = terms.agrid;
            Kgrid = terms.Kgrid;
            pil = terms.pil;
            g = terms.G;
            captax = terms.captax;
            Kfore = terms.Kfore;
            p = terms.p;

            VP = zeros(nm, nl, na);
            VL = VP;
            GP = VP;
            GL = GP;
            TVP = VP;
            TVL = VL;
            TGP = GP;
            TGL = GL;

            V0 = zeros(nl, na);

            % set up V so that it doesn't start empty
            % using believable r and w values: 4% interest, w = 1.3
            scale = .25;
            for im = 1:nm
                for ia = 1:na
                    kval = agrid(ia);
                    for il = 1:nl
                        yval = scale*(1+0.04*(1-captax(il)))*kval + ...
                            1.3*lgrid(il) - 0.04*phi;
                        ymin = max(1e-10, yval);
                        VL(im, il, ia) = log(ymin);
                    end
                end
            end
            VP = VL;

            %init expected vals

            V = VP;
            EV = ks.getExpectation(V,Kgrid, Kfore, p);
            for im = 1:nm
                for ia = 1:na
                    for il = 1:nl
                        V0(im, il, ia) = pil(il,:)*EV(im, :,ia)';
                    end
                end
            end

            dist = 1e5;
            iter_ct = 1;
            
            while dist > vTol
                
                EVP = ks.getExpectation(VP,Kgrid, Kfore, p);
                EVL = ks.getExpectation(VL,Kgrid, Kfore, p);
                for im = 1:nm
                    for ia = 1:na
                        for il = 1:nl
                            V0(im, il, ia) = p*pil(il,:)*EVP(im, :,ia)' ...
                                + (1-p)*pil(il,:)*EVL(im, :,ia)';
                        end
                    end
                end

                % need to forecast prices and then put them into EGM for
                % each moment. Then need to store the expectation and use
                % that to calculate ...? this is in the getregdata part
                CP = TVP; CL = TVL;
                for i = 1:nm
                    K = Kgrid(i);
                    Kp = ks.forecast(Kfore(1,:), K);
                    Kl = ks.forecast(Kfore(2,:), K);
                    terms.r = vaas.calcr(terms.alpha, terms.delta, Kp, terms.etap);
                    terms.w = vaas.calcw(terms.alpha, Kp, terms.etap);
                    [TVP(i,:,:) TGP(i,:,:) CP(i,:,:)]= egm.solve(nl, na, terms,squeeze(V0(i,:,:)), ...
                        squeeze(VP(i,:,:)));
                    terms.r = vaas.calcr(terms.alpha, terms.delta, Kl, terms.etal);
                    terms.w = vaas.calcw(terms.alpha, Kl, terms.etal); 
                    [TVL(i,:,:) TGL(i,:,:) CL(i,:,:)]= egm.solve(nl, na, terms,squeeze(V0(i,:,:)), ...
                        squeeze(VL(i,:,:)));
                end 

                dist = max(compute.dist(VP, TVP, 3), compute.dist(VL, TVL, 3));
                kdist = max(compute.dist(GP, TGP, 3), compute.dist(GL, TGL, 3));
            
                if mod(iter_ct, 20) == 0
                    fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
                        "\n\t\t||TG - G|| = %4.6f", iter_ct, dist, kdist);
                    fprintf("\nInitial Values:");
                    fprintf("\nMin Kgrid: %2.4f, Max Kgrid: %2.4f", min(Kgrid), max(Kgrid));
                    fprintf("\nInitial Capital Forecasts: Kp = %2.4f, Kl = %2.4f", Kp, Kl);
                    fprintf("\nInitial Policy Function: Min Gp = %2.4f, Max Gp = %2.4f", ...
                        min(TGP(:)), max(TGP(:)));
                     fprintf("\nInitial Policy Function: Min Gl = %2.4f, Max Gl = %2.4f", ...
                        min(TGL(:)), max(TGL(:)));
                   fprintf("\nInitial Value Function: Min VP = %2.4f, Max VP = %2.4f", ...
                       min(TVP(:)), max(TVP(:)));
                   fprintf("\nInitial Value Function: Min VL = %2.4f, Max VL = %2.4f", ...
                       min(TVL(:)), max(TVL(:)));
                end
            
                iter_ct = iter_ct + 1;
            
                VP = TVP;
                VL = TVL;

            end

            if verbose
                fprintf("\n\tIteration %i: ||TV - V|| = %4.6f" + ...
                    "\t||TG - G|| = %4.6f\n", iter_ct, dist, kdist);
            end

        end

        %% get expectation
        % V is a three-dimensional object, because it's also defined on the
        % aggregate state grid. This function approximates expecations of V
        % given the aggregate state grid by linearly interpolating the
        % expecation over the aggregate states. The forecasting rule gives
        % our expected future state.
        function EV = getExpectation(V, Kgrid, Kfore, p)

            nm = length(Kgrid);

            EV = zeros(size(V));
            for i = 1:nm
                Ktoday = Kgrid(i);
                K1 = ks.forecast(Kfore(1,:), Ktoday);
                K2 = ks.forecast(Kfore(2,:), Ktoday);
                Ktomorrow = p*K1 + (1-p)*K2;
                [ix, we] = compute.weight(Kgrid, Ktomorrow);
                EV(i, :, :) = we*V(ix, :, :) + (1-we)*V(ix+1, :, :);
            end

        end

        function kpr = forecast(fore,k)

            b0 = fore(1,1);
            b1 = fore(1,2);
            kpr = exp(b0)*(k^b1);
        end

    end
end