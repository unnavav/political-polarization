classdef HH
    methods(Static)

        %% solve the HH problem
        % this function is used four times to solve the four different 
        % HH x Party types : (a,A), (b,A), (a,B), (b,B). The things that
        % change between combos:
        %   - ubonus
        %   - tau
        %   - g
        function [V, G, C, V0] = solve(nl, na, terms, vTol, verbose)

            beta = terms.beta;
            sigma = terms.sigma;
            phi = terms.phi;
            lgrid = terms.lgrid;
            agrid = terms.agrid;
            pil = terms.pil;
            g = terms.G;
            captax = terms.captax;
            lamval = terms.lamval;
            tau = terms.tau;
            p = terms.p;

            r = terms.r;
            w = terms.w;

            V = zeros(nl, na);
            G = V;
            TG = G;
            TV = V; 

            V0 = zeros(nl, na);

            % set up V so that it doesn't start empty
            scale = .25;
            for ia = 1:na
                kval = agrid(ia);
                for il = 1:nl
                    yval = scale*(1+r*(1-captax(il)))*kval + w*lgrid(il) - r*phi;
                    ymin = max(1e-10, yval);
                    V(il, ia) = log(ymin);
                end
            end
            
            %init expected vals
            for ia = 1:na
                for il = 1:nl
                    V0(il, ia) = pil(il,:)*V(:,ia);
                end
            end

            dist = 1e5;
            iter_ct = 1;
            
            while dist > vTol
                        
                %endogrid choice
                endoK = zeros(size(V));
                for ia = 1:na
                    kpr = agrid(ia);
                    for il = 1:nl
                        l = lgrid(il);
                        D = egm.solveD(V0(il,:), ia, agrid);

                        endoK(il, ia) = ((beta*D)^(-1) + kpr - w*l +...
                            gov.tax(w*l, lamval, tau) - g(1) + ...
                            r*(1-captax(il))*phi)/(1+r*(1-captax(il)));
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
            
                C = endoK;
                for ia = 1:na
                    for il = 1:nl
                        l = lgrid(il);

                        c = (1+r*(1-captax(il)))*agrid(ia) + w*l - gov.tax(w*l, lamval, tau) ...
                            + g - G(il, ia) - r*(1-captax(il))*phi;

                        C(il, ia) = max(1e-10, c);

                        [ix, we] = compute.weight(agrid, G(il, ia));
                        ev = we*V0(il, ix) + (1-we)*V0(il, ix+1);
                        TV(il, ia) = log(c) + beta*ev;
                    end
                end
            
                dist = compute.dist(TV, V, 2);
                kdist = compute.dist(TG, G, 2);
            
%                 if mod(iter_ct, 50) == 0
%                     fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
%                         "\n\t\t||TG - G|| = %4.6f", iter_ct, dist, kdist);
%                 end
            
                iter_ct = iter_ct + 1;
            
                V = TV;
                TG = G;

                for ia = 1:na
                    for il = 1:nl
                        V0(il, ia) = pil(il,:)*V(:,ia);
                    end
                end

            end

            if verbose
                fprintf("\n\tIteration %i: ||TV - V|| = %4.6f" + ...
                    "\t||TG - G|| = %4.6f\n", iter_ct, dist, kdist);
            end

        end
        
        % the issue is that I can't use EGM to solve, I need to use GSS.
        % sad reacts. So implementing that here, given that E0 is already
        % set, so I can get smooth policy functions.
        function [V, G] = backsolve(nl, na, Vpr, terms, vTol, verbose)

            beta = terms.beta;
            sigma = terms.sigma;
            phi = terms.phi;
            lgrid = terms.lgrid;
            agrid = terms.agrid;
            pil = terms.pil;
            g = terms.G;
            captax = terms.captax;
            lamval = terms.lamval;
            tau = terms.tau;

            r = terms.r;
            w = terms.w;

            V = zeros(nl, na);
            G = V;
            TG = G;
            TV = V; 

            % set up V so that it doesn't start empty
            scale = .25;
            for ia = 1:na
                kval = agrid(ia);
                for il = 1:nl
                    yval = scale*(1+r*(1-captax(il)))*kval + w*lgrid(il) - r*phi;
                    ymin = max(1e-10, yval);
                    V(il, ia) = log(ymin);
                end
            end
            
            % expected vals
            for ia = 1:na
                for il = 1:nl
                    V0(il, ia) = pil(il,:)*Vpr(:,ia);
                end
            end
            
            % get best capital choice from V, EV using GSS
            for ia = 1:na

                for il = 1:nl
                    y = w*lgrid(il) + (1+r*agrid(ia));
                    Cvals = V0(il, :);
                    params = [y beta sigma];

                    [v, g, ~] = compute.gss(Cvals, params, agrid, vTol/1e3);

                    V(il, ia) = v;
                    G(il, ia) = g;
                end
            end
        
            dist = compute.dist(TV, V, 2);
            kdist = compute.dist(TG, G, 2);
        
            if verbose
                fprintf("\n\tIteration %i: ||TV - V|| = %4.6f" + ...
                    "\t||TG - G|| = %4.6f\n", iter_ct, dist, kdist);
            end
        end

        % we have two types of heterogeniety: wealth and political agent
            % must solve for distributions of each type by which party is in
            % power, and then scale the distribution by the percentage of
            % HHs at that point who are of type A
        % e.g. 25% of HH a's have savings of 100 under party A, but
        % only 50% of HH are party a, so in actuality only 12.5% of
        % distr is at 100. 
        function [mu1, kagg] = getDist(G, amu, agrid, pil, verbose)

            [nl, ~] = size(G);
            [~, nmu] = size(amu);
            mu = zeros(nl, nmu);

            %nomenclature: ixagrid is indices for HH a for both parties.
            ixgrid = zeros(size(G)); wegrid = ixgrid;
          
            % We linearly interpolate the policy function on agrid to compute
            % afval on the distribution support.
            for im = 1:nmu
                kval = amu(im);
                for il = 1:nl  
                    [ix, we] = compute.weight(agrid, kval);
                        
                    %split between rep and dem capital choices
                    kdval = G(il,ix)*we + G(il,ix+1)*(1.0 - we);

                    [ix, we] = compute.weight(amu, kdval);
                    ixgrid(il, im) = ix;
                    wegrid(il, im) = we;
                end
            end
            
            distance = 20; iter_ct = 1;
            
            % nomenclature: muA is the distribution over assets when A is 
            % in power
            mu = ones(size(mu))*(1/(nmu*nl));

            while (distance > 1e-6 && iter_ct < 3000)
                
                mu1 = zeros(size(mu));
                
                for im = 1:nmu
                    for il = 1:nl
                        
                        %need weighting for a and b HH under party A
                        ix = ixgrid(il,im); we = wegrid(il,im); 
                        muval = mu(il, im); 
                        
                        if muval > 0
                            for jl = 1:nl

                                %a households' movements
                                if (ix < nmu)
                                    mu_val1 = pil(il,jl)*muval*we;
                                    mu_val2 = pil(il,jl)*muval*(1.0 - we);
                                else
                                    mu_val1 = pil(il,jl)*muval*we;
                                    mu_val2 = 0;
                                end

                                mu1(jl,ix) = mu1(jl,ix) + mu_val1;
                                mu1(jl,ix+1) = mu1(jl,ix+1) + mu_val2;
                            end
                        end
                    end
                end
                                
                distance = compute.dist(mu1, mu, 2);

%                 if (mod(iter_ct,50) == 0)
%                     s = sprintf( '\n\t\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
%                         iter_ct, distance, sum(sum(mu1)));
%                     disp(s)
%                 end
                
                iter_ct = iter_ct + 1;
            
                mu = mu1;   
            end
            s = sprintf( '\n\t\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
                iter_ct, distance, sum(sum(mu)));
            if verbose
                disp(s);
            end 
            
            distrA2500 = sum(mu,1);
            kagg = amu*distrA2500';
        end

        function [mu1, kagg] = transitDistr(G, mu, amu, agrid, pil)

            [nl, ~] = size(G);
            [~, nmu] = size(amu);

            %nomenclature: ixagrid is indices for HH a for both parties.
            ixgrid = zeros(size(G)); wegrid = ixgrid;
          
            % We linearly interpolate the policy function on agrid to compute
            % afval on the distribution support.
            for im = 1:nmu
                kval = amu(im);
                for il = 1:nl  
                    [ix, we] = compute.weight(agrid, kval);
                        
                    %split between rep and dem capital choices
                    kdval = G(il,ix)*we + G(il,ix+1)*(1.0 - we);

                    [ix, we] = compute.weight(amu, kdval);
                    ixgrid(il, im) = ix;
                    wegrid(il, im) = we;
                end
            end
                
            mu1 = zeros(size(mu));
            
            for im = 1:nmu
                for il = 1:nl
                    
                    ix = ixgrid(il,im); we = wegrid(il,im); 
                    muval = mu(il, im); 
                    
                    if muval > 0
                        for jl = 1:nl

                            %a households' movements
                            if (ix < nmu)
                                mu_val1 = pil(il,jl)*muval*we;
                                mu_val2 = pil(il,jl)*muval*(1.0 - we);
                            else
                                mu_val1 = pil(il,jl)*muval*we;
                                mu_val2 = 0;
                            end

                            mu1(jl,ix) = mu1(jl,ix) + mu_val1;
                            mu1(jl,ix+1) = mu1(jl,ix+1) + mu_val2;
                        end
                    end
                end
            end
                              
            distrA2500 = sum(mu,1);
            kagg = amu*distrA2500';
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

            if sigma == 1
                utils = log(y-apr) + beta*vc;
            else
                utils = ((y-apr)^(1-sigma))/(1-sigma) + beta*vc;
            end
        end
    end
end