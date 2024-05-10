classdef HH
    methods(Static)


        function [V, g] = solve(nl, na, terms, vTol)

            beta = terms.beta;
            sigma = terms.sigma;
            phi = terms.phi;
            r = terms.r;
            wage = terms.wage;
            lgrid = terms.lgrid;
            agrid = terms.agrid;
            pil = terms.pil;
            lambda = terms.lambda;
            tau = terms.tau;

            V = zeros(nl, na);
            EV = V;
            g = V;
            V1 = V;
            E = EV;

            % initialize value functions
            V = zeros(nl, na); %(labor, assets, my type, goverment type)

            for il = 1:nl
                l = lgrid(il);
                for ia = 1:na
                    a = agrid(ia);
            
                    max_cons = wage*l + (1+r)*a - r*phi;
                       
                    max_cons = max_cons - T(il, ia);
                         
                    if sigma == 1
                        V(il, ia) = log(max_cons)/(1-beta); 
                    else
                        V(il, ia) = (max_cons^(1-sigma))/((1-beta)*(1-sigma));
                    end

                end
            end

            dist = 10;
            iter_ct = 1;
            options = optimoptions('fsolve','Display','none');
            
            while dist > vTol
            
                % get expected value function to start with
                for il = 1:nl
                    EV(il, :) = pil(il,:)*V(:,:);
                end

                DEV = egm.numdev(EV,agrid);

                %now find endogenous grid
            
                parfor il = 1:nl
                    w_e = wage*lgrid(il);
                    for apr = 1:na            
                        a_pr = agrid(apr);
                        c = (beta*DEV(il, apr))^(-1/sigma);
                        capr = c+w_e;
                        
                        bcterms = struct('we', w_e, ...
                            'capr', capr, ...
                            'r', r, ...
                            'tau', tau);
                        budgetfun = @(x)egm.impliedbc(x, bcterms);
                        E(il, ia)  = fsolve(budgetfun, a_pr, options);      %%%% see if this slows it down        
                        

                    end
                end
                
                % now project it back onto agrid
    
                for il = 1:nl
                    w_e = wage*lgrid(il);
                    
                    for ia = 1:na
                        lb = min(E(il, :));
                        
                        ahat = agrid(ia);
                        
                        if ahat < lb
                            g(il, ia) = 0;
                        else
                            [ix, we] = compute.weight(E(il, :), ahat);
                            g(il, ia) = we*E(il, ix) + (1-we)*E(il, ix + 1);
                        end
                        
                        c = gov.tax((1+r)*ahat + w_e, lambda, tau) ...
                            - r*phi - g(il, ia);
                        
                        if c < 0
                            disp([il ia])
                        end
                        
                        if sigma == 1
                            V(il, ia) = log(c) + beta*V1(il, ia);
                        else
                            V(il, ia) = (c^(1-sigma))/(1-sigma) + ...
                                beta*V1(il, ia);
                        end
                        
                    end
                end
            end

            dist = compute.dist(V, V1, 2);
        
%                 if mod(iter_ct, 100) == 0
%                     fprintf("\n\t\tIteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);
%                 end
        
            iter_ct = iter_ct + 1;
        
            V1 = V;            
        end

        fprintf("\n\t\tIteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);

    end

        function [mu1, Kagg] = getDist(g, amu, agrid, pil)

            [nl, ~] = size(g);
            nmu = max(size(amu));
            mu = zeros(nl, nmu); ixgrid = mu; wegrid = mu;

            for im = 1:nmu
                kval = amu(im);
                for il = 1:nl  
                        % We linearly interpolate the policy function on bgrid to compute
                        % bfval on the distribution support.
                        
                        [ix, we] = compute.weight(agrid, kval);
                            
                        %split between rep and dem capital choices
                        kdval = g(il,ix)*we + g(il,ix+1)*(1.0 - we);
                        
                        [ix, we] = compute.weight(amu, kdval);
                                
                        ixgrid(il, im) = ix;
                        wegrid(il, im) = we;
                end
            end
            
            distance = 20; iter_ct = 1;
            
            for i = 1:nmu
                for il = 1:nl
                    mu(il,i) = 1.0/(nmu*nl);
                end
            end

            while (distance > 1e-8 && iter_ct < 3000)
                
                mu1 = zeros(size(mu));
                
                for im = 1:nmu
                    for il = 1:nl

                        ix = ixgrid(il,im); we = wegrid(il,im); muval = mu(il,im);
                        
                        if muval > 0
                            for jl = 1:nl
                                if (ix < nmu)
                                    mu1(jl,ix) = mu1(jl,ix) + pil(il,jl)*muval*we;
                                    mu1(jl,ix+1) = mu1(jl,ix+1) + pil(il,jl)*muval*(1.0 - we);
                                else
                                    mu1(jl,ix) = mu1(jl,ix) + pil(il,jl)*muval*we;
                                end
                            end
                        end

                    end
                end
                                
                distance = compute.dist(mu1, mu, 2);
                
%                 if (mod(iter_ct,200) == 0)
%                     s = sprintf( '\n\t\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
%                         iter_ct, distance, sum(sum(mu1)));
%                     disp(s);
%                 end
                
                iter_ct = iter_ct + 1;
            
                mu = mu1;   
            end
            s = sprintf( '\n\t\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
                iter_ct, distance, sum(sum(mu)));
            disp(s);

            distr2500 = sum(mu,1);
            Kagg = amu*distr2500';

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

    end
end