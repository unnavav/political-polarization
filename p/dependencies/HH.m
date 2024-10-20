classdef HH
    methods(Static)

        %% solve the HH problem
        % this function is used four times to solve the four different 
        % HH x Party types : (a,A), (b,A), (a,B), (b,B). The things that
        % change between combos:
        %   - ubonus
        %   - tau
        %   - g
        function [V, G, V0] = solve(nl, na, np, terms, vTol)

            beta = terms.beta;
            sigma = terms.sigma;
            phi = terms.phi;
            K = terms.K;
            L = terms.L;
            lgrid = terms.lgrid;
            agrid = terms.agrid;
            pil = terms.pil;
            g = terms.G;
            ubonus = terms.ubonus;
            captax = terms.captax;
            lamval = terms.lamval;
            tau = terms.tau;
            p = terms.p;
            eta = terms.eta;

            % no migration (populist)
            rp = alpha*(K^(alpha - 1)*(L^(1-alpha))) - delta;
            wp = (1-alpha)*((K^(alpha))*(L^(-alpha)));

            % yes migration (liberalist)
            rl = alpha*(K^(alpha - 1)*((L*eta)^(1-alpha))) - delta;
            wl = (1-alpha)*((K^(alpha))*((L*eta)^(-alpha)));

            VL = zeros(nl, na);
            VP = zeros(nl, na);
            G = VL;
            TV = VL; V = VL; 

            V0 = zeros(nl, na);

            % set up V so that it doesn't start empty
            scale = .0025;
            for ia = 1:na
                kval = agrid(ia);
                for il = 1:nl
                    yval = scale*(1+r*(1-captax(il)))*kval + w*lgrid(il) - r*phi;
                    ymin = max(1e-10, yval);
                    VL(il, ia) = log(ymin);
                    VP(il, ia) = log(ymin);
                end
            end
            
            %init expected vals
            for ia = 1:na
                for il = 1:nl
                    for ip = 1:np
                        prob = p(ip);
                        V0(il, ia, ip) = prob*pil(il,:)*VP(:,ia) + ...
                            (1-prob)*pil(il,:)*VB(:,ia);
                    end
                end
            end

            dist = 1e5;
            iter_ct = 1;
            
            while dist > vTol
            
                % making a change--I calculate c_{t+1} given our decision
                % rule at productivity il, asset ia. then use that to get
                % expected value. So will need a different VA,VB   
                for ia = 1:na
                    for il = 1:nl
                        ca = (1+r)*agrid(ia) + gov.tax(wl*lgrid(il), lamval, ...
                            tau(1)) - G(il, ia);
                        VP(il, ia) = log(ca) + beta*V0(il, ia);
                        cb = (1+r2)*agrid(ia) + gov.tax(wp*lgrid(il), lamval, ...
                            tau(2)) - G(il, ia);
                        VL(il, ia) = log(cb) + beta*V0(il, ia);
                    end
                end
            
                %now get expectation
                for ia = 1:na
                    for il = 1:nl
                        V0(il, ia) = p*pil(il,:)*VP(:,ia) + ...
                            (1-prob)*p(il,:)*VB(:,ia);
                    end
                end

                %endogrid choice: back out K from consumption of EE
                endoK = zeros(size(VP));
                for ia = 1:na
                    kpr = agrid(ia);
                    for il = 1:nl
                        l = lgrid(il);                        
                        c = exp(beta*(p*VP(il, ia)+(1-p)*VL(il,ia)));
                        endoK(il, ia) = (1+rl)*kpr + gov.tax(wl*l,lamval, ...
                            tgrid(1+round(p/1)))-c;
                    end
                end
            
                %linterpolate decision rule
                for ia = 1:na
                    k = agrid(ia);
                    for il = 1:nl
                        lkvals = endoK(il, :);
                        
                        if k < lkvals(1)
                            G(il, ia) = agrid(1);
                        else 
                            [ix, we] = compute.weight(lkvals, k);
                            kpr = we*agrid(ix) + (1-we)*agrid(ix + 1);
                            G(il, ia) = max(agrid(1), kpr);
                        end

                    end
                end
            
                for ia = 1:na
                    for il = 1:nl
                        l = lgrid(il);

                        %A in power
                        c = (1+r*(1-captax(il)))*agrid(ia) + w*l - gov.tax(w*l, lamval, tau(1)) ...
                            + g(1) - G(il, ia) - r*(1-captax(il))*phi;
                        [ix, we] = compute.weight(agrid, G(il, ia));
                        ev = we*V0(il, ix,1) + (1-we)*V0(il, ix+1,1);
                        V(il, ia) = log(c) + beta*ev;

                    end
                end
            
                dist = compute.dist(TV, V, 2);
                kdist = max(compute.dist(TG, G, 2));
            
%                 if mod(iter_ct, 10) == 0
%                     fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
%                         "\n\t\t||TG - G|| = %4.6f", iter_ct, dist, kdist);
%                 end
            
                iter_ct = iter_ct + 1;
            
                TV = V;
                TG = G;
            end

        
            fprintf("\n\tIteration %i: ||TV - V|| = %4.6f" + ...
                "\t||TG - G|| = %4.6f\n", iter_ct, dist, kdist);

            V(:,:,1) = VA; V(:,:,2) = VB;
            G(:,:,1) = GA; G(:,:,2) = GB;
        end

        % we have two types of heterogeniety: wealth and political agent
            % must solve for distributions of each type by which party is in
            % power, and then scale the distribution by the percentage of
            % HHs at that point who are of type A
        % e.g. 25% of HH a's have savings of 100 under party A, but
        % only 50% of HH are party a, so in actuality only 12.5% of
        % distr is at 100. 
        function [mu1A, KaggA, mu1B, KaggB] = getDist(ga, gb, amu, agrid, pil, pctA)

            [nl, ~, np] = size(ga);
            [~, nmu] = size(amu);
            muA = zeros(nl, nmu);
            muB = zeros(nl, nmu);

            %nomenclature: ixagrid is indices for HH a for both parties.
            ixagrid = zeros(size(ga)); weagrid = ixagrid; 
            ixbgrid = ixagrid; webgrid = ixagrid;
          
            % We linearly interpolate the policy function on agrid to compute
            % afval on the distribution support.
            for im = 1:nmu
                kval = amu(im);
                for il = 1:nl  
                    for ip = 1:np

                        
                        [ix, we] = compute.weight(agrid, kval);
                            
                        %split between rep and dem capital choices
                        kdaval = ga(il,ix, ip)*we + ga(il,ix+1, ip)*(1.0 - we);
                        kdbval = gb(il,ix, ip)*we + gb(il,ix+1, ip)*(1.0 - we);

                        [ixa, wea] = compute.weight(amu, kdaval);
                        ixagrid(il, im, ip) = ixa;
                        weagrid(il, im, ip) = wea;

                        [ixb, web] = compute.weight(amu, kdbval);
                        ixbgrid(il, im, ip) = ixb;
                        webgrid(il, im, ip) = web;
                    end
                end
            end
            
            distance = 20; iter_ct = 1;
            
            % nomenclature: muA is the distribution over assets when A is 
            % in power
            for im = 1:nmu
                for il = 1:nl
                    muA(il,im) = 1.0/(nmu*nl);
                    muB(il,im) = 1.0/(nmu*nl);
                end
            end

            while (distance > 1e-6 && iter_ct < 3000)
                
                mu1A = zeros(size(muA));
                mu1B = zeros(size(muB));
                
                for im = 1:nmu
                    for il = 1:nl
                        
                        %Party A ---------------
                        ip = 1;

                        %need weighting for a and b HH under party A
                        ix_a = ixagrid(il,im, ip); we_a = weagrid(il,im, ip); 
                        ix_b = ixbgrid(il,im, ip); we_b = webgrid(il,im, ip); 
                        muvalA = muA(il, im); 
                        
                        if muvalA > 0
                            for jl = 1:nl

                                %a households' movements
                                if (ix_a < nmu)
                                    hha_mu1 = pil(il,jl)*muvalA*we_a;
                                    hha_mu2 = pil(il,jl)*muvalA*(1.0 - we_a);
                                else
                                    hha_mu1 = pil(il,jl)*muvalA*we_a;
                                    hha_mu2 = 0;
                                end

                                %b households' movements
                                if (ix_b < nmu)
                                    hhb_mu1  = pil(il,jl)*muvalA*we_b;
                                    hhb_mu2 = pil(il,jl)*muvalA*(1.0 - we_b);
                                else
                                    hhb_mu1 = pil(il,jl)*muvalA;
                                    hhb_mu2 = 0;
                                end

                                mu1A(jl,ix_a) = mu1A(jl,ix_a) + pctA*hha_mu1 ...
                                    + (1-pctA)*hhb_mu1;
                                mu1A(jl,ix_a+1) = mu1A(jl,ix_a+1) + pctA*hha_mu2 ...
                                    + (1-pctA)*hhb_mu2;
                            end
                        end

                        %Party B ---------------
                        ip = 2;
                        ix_a = ixagrid(il,im, ip); we_a = weagrid(il,im, ip); 
                        ix_b = ixbgrid(il,im, ip); we_b = webgrid(il,im, ip);
                        muvalB = muB(il, im);

                        if muvalB > 0
                            for jl = 1:nl

                                %a households' movements
                                if (ix_a < nmu)
                                    hha_mu1 = pil(il,jl)*muvalB*we_a;
                                    hha_mu2 = pil(il,jl)*muvalB*(1.0 - we_a);
                                else
                                    hha_mu1 = pil(il,jl)*muvalB*we_a;
                                    hha_mu2 = 0;
                                end

                                %b households' movements
                                if (ix_b < nmu)
                                    hhb_mu1  = pil(il,jl)*muvalB*we_b;
                                    hhb_mu2 = pil(il,jl)*muvalB*(1.0 - we_b);
                                else
                                    hhb_mu1 = pil(il,jl)*muvalB;
                                    hhb_mu2 = 0;
                                end

                                mu1B(jl,ix_a) = mu1B(jl,ix_a) + pctA*hha_mu1 ...
                                    + (1-pctA)*hhb_mu1;
                                mu1B(jl,ix_a+1) = mu1B(jl,ix_a+1) + pctA*hha_mu2 ...
                                    + (1-pctA)*hhb_mu2;
                            end
                        end

                    end
                end
                                
                distanceA = compute.dist(mu1A, muA, 3);
                distanceB = compute.dist(mu1B, muB, 3);

%                 if (mod(iter_ct,50) == 0)
%                     s = sprintf( '\n\t\tIteration %3i: ||Tma - ma|| = %8.6f\tsum = %6.4f ', ...
%                         iter_ct, distanceA, sum(sum(mu1A)));
%                     disp(s);
%                     s = sprintf( '\n\t\t\t\t\t\t||Tmb - mb|| = %8.6f\tsum = %6.4f ', ...
%                         distanceB, sum(sum(mu1B)));
%                     disp(s);
%                 end
                
                distance = max(distanceA, distanceB);

                iter_ct = iter_ct + 1;
            
                muA = mu1A;   
                muB = mu1B;   
            end
            s = sprintf( '\n\t\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
                iter_ct, distance, sum(sum(muA)));
            disp(s);

            distrA2500 = sum(muA,1);
            KaggA = amu*distrA2500';

            distrB2500 = sum(muB,1);
            KaggB = amu*distrB2500';
        end


        function [vdistr, winner] = map(VOTES, amu, agrid, adistr, pctDem)

            [nl, np, nm]=size(adistr);
            ixgrid = zeros(nm); wegrid = ixgrid;
            vdistr = zeros(nl, nm);
            %first get mapping
            for im=1:nm
                a = amu(im);
                [ix, we] = compute.weight(agrid, a);

                ixgrid(im) = ix; wegrid(im) = we;
            end

            for im = 1:nm
                for ip = 1:np
                    for il = 1:nl

                        muval = adistr(il, ip, im);
                        if muval>0
                            ix = ixgrid(im); we = wegrid(im);
                            vdistr(il, im) = pctDem*VOTES(il, ix,1)*we*muval+...
                                (1-pctDem)*VOTES(il, ix,2)*we*muval + ...
                                pctDem*VOTES(il, ix+1,1)*(1-we)*muval+...
                                (1-pctDem)*VOTES(il, ix+1,2)*(1-we)*muval;
                        end
                    end
                end
            end

            majority = sum(sum(vdistr));

            if majority > .5
                winner = 1;
            else 
                winner = 0;
            end
        end

        function utils = util(apr, params, vc)
            y = params(1);
            beta = params(2);
            sigma = params(3);
            ubonus = params(4);

            if y-apr <= 0
                utils = -100000;
            elseif sigma == 1
                utils = log(y-apr) + ubonus + beta*vc;
            else
                utils = ((y-apr)^(1-sigma))/(1-sigma) + ubonus + beta*vc;
            end
        end
    end
end