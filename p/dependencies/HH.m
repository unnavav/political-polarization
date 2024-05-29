classdef HH
    methods(Static)


        function [V, G, V0] = solve(nl, na, np, terms, vTol, gTol)

            beta = terms.beta;
            sigma = terms.sigma;
            phi = terms.phi;
            r = terms.r;
            w = terms.wage;
            lgrid = terms.lgrid;
            agrid = terms.agrid;
            pil = terms.pil;
            g = terms.G;
            p = terms.p;
            ubonus = terms.ubonus;
            captax = terms.captax;
            lamval = terms.lamval;
            tau = terms.tau;
            
            V = zeros(nl, na);
            V0 = V;
            TV = V;
            
            G = V;
            TG = V;
            
            % set up V so that it doesn't start empty
            scale = 0.025;
            for ia = 1:na
                kval = agrid(ia);
                for il = 1:nl
                    yval = scale*(1+r*(1-captax))*kval + w*lgrid(il) - r*phi;
                    ymin = max(1e-10, yval);
                    V(il, ia) = log(ymin);
                end
            end
            
            
            dist = 1e5;
            iter_ct = 1;
            
            while dist > vTol
            
                for ia = 1:na
                    for il = 1:nl
                        V0(il, ia) = pil(il,:)*V(:,ia);
                    end
                end
            
                %endogrid choice
                endoK = zeros(size(V));
                for ia = 1:na
                    kpr = agrid(ia);
                    for il = 1:il
                        l = lgrid(il);
                        D = egm.solveD(V0(il,:), ia, agrid);
                        endoK(il, ia) = ((beta*D)^(-1) + kpr - ...
                            gov.tax(w*l, lamval, tau) - g + ...
                            r*(1-captax)*phi)/(1+r*(1-captax));
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
                        c = (1+r*(1-captax))*agrid(ia) + gov.tax(w*l, lamval, tau) ...
                            + g - G(il, ia) - r*(1-captax)*phi;
                        [ix, we] = compute.weight(agrid, G(il, ia));
                        ev = we*V0(il, ix) + (1-we)*V0(il, ix+1);
                        TV(il, ia) = log(c) + beta*ev;
                    end
                end
            
                dist = compute.dist(TV, V, 2);
                kdist = compute.dist(TG, G, 2);
            
                if mod(iter_ct, 200) == 0
                    fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
                        "\n\t\t||TG - G|| = %4.6f", iter_ct, dist, kdist);
                end
            
                iter_ct = iter_ct + 1;
            
                V = TV;
                TG = G;
            end

        
            fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
                "\n\t\t||TG - G|| = %4.6f\n", iter_ct, dist, kdist);

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

                        [ix, we] = compute.weight(amu, kdaval);
                        ixagrid(il, im, ip) = ix;
                        weagrid(il, im, ip) = we;

                        [ix, we] = compute.weight(amu, kdbval);
                        ixbgrid(il, im, ip) = ix;
                        webgrid(il, im, ip) = we;
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

                if (mod(iter_ct,50) == 0)
                    s = sprintf( '\n\t\tIteration %3i: ||Tma - ma|| = %8.6f\tsum = %6.4f ', ...
                        iter_ct, distanceA, sum(sum(mu1A)));
                    disp(s);
                    s = sprintf( '\n\t\t\t\t\t\t||Tmb - mb|| = %8.6f\tsum = %6.4f ', ...
                        distanceB, sum(sum(mu1B)));
                    disp(s);
                end
                
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