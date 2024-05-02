classdef HH
    methods(Static)


        function [V, g, VOTES] = solve(nl, na, np, terms, vTol)

            beta = terms.beta;
            sigma = terms.sigma;
            phi = terms.phi;
            identity = terms.identity;
            r = terms.r;
            wage = terms.wage;
            lgrid = terms.lgrid;
            agrid = terms.agrid;
            tgrid = terms.tgrid;
            lamgrid = terms.lamgrid;
            pil = terms.pil;
            p = terms.p;

            V = zeros(nl, na, np, np);
            EV = V;
            g = V;
            V1 = V;

            VOTES = zeros(nl, na, np);

            % initialize value functions
            V = zeros(nl, na, np, np); %(labor, assets, my type, goverment type)
            
            for il = 1:nl
                l = lgrid(il);
                for ia = 1:na
                    a = agrid(ia);
            
                    max_cons = wage*l + (1+r)*a - r*phi;
                        
                    for ip = 1:np
            
                        max_cons = gov.tax(max_cons, lamgrid(ip), tgrid(ip));
                        
                        if sigma == 1
                            V(il, ia, ip, ip) = log(max_cons)/(1-beta); 
                        else
                            V(il, ia, ip, ip) = (max_cons^(1-sigma))/((1-beta)*(1-sigma));
                        end
                    end
                end
            end

            % get expected value function to start with
            for il = 1:nl
                for ip = 1:np
                    EV(il, :, 1, ip) = (p*pil(il,:)*V(:,:,1,1) + ...
                        (1-p)*pil(il,:)*V(:,:,1,2))';
                    EV(il, :, 2, ip) = (p*pil(il,:)*V(:,:,2,1) + ...
                        (1-p)*pil(il,:)*V(:,:,2,2))';
                end
            end

            dist = 10;
            iter_ct = 1;
            
            while dist > vTol
            
                %g choice: maximize expected value
            
                for il = 1:nl
                    l = lgrid(il);
                    for ia = 1:na
            
                        a = agrid(ia);
            
                        y = (1+r)*a + wage*l - r*phi;
            
                        yvec = ones(250,1)*y;
            
                        C = yvec - agrid';
                        C = C(C>0);
                        C2 = C;
                
                        for ip = 1:np
            
                            nchoices = max(size(C));
                            for ic = 1:nchoices
                                C2(ic) = gov.tax(C(ic),lamgrid(ip), tgrid(ip));
                            end
            
                            %C = (1-tgrid(ip))*C;
            
                            C = C2;

                            if sigma == 1
                                Val = log(C)' + beta*EV(il, 1:size(C,1),ip);
                            else
                                Val = (C.^(1-sigma))/(1-sigma)' + beta*EV(il, 1:size(C,1), ip);
                            end
                            
                            [val, ai] = max(Val);
                
                            if ip - 1 > 0
                                V(il, ia, 1, ip) = val;
                                g(il, ia, 1, ip) = agrid(ai);
                                V(il, ia, 2, ip) = val + identity;
                                g(il, ia, 2, ip) = agrid(ai);
                            else 
                                V(il, ia, 1, ip) = val + identity;
                                g(il, ia, 1, ip) = agrid(ai);
                                V(il, ia, 2, ip) = val;
                                g(il, ia, 2, ip) = agrid(ai);
                            end
                        end
                    end
                end
            
                % get expected value function after getting value function
            
                %TODO ADD PARTY DIFFERENCES
                for il = 1:nl
                    for ip = 1:np
                        EV(il, :, 1, ip) = pil(il,:)*V(:,:,1,ip);
                        EV(il, :, 2, ip) = pil(il,:)*V(:,:,2,ip);
                    end
                end
            
                % calculate voting decision
            
                VOTES(:,:,1) = EV(:,:,1,1) >= EV(:,:,1,2);
                VOTES(:,:,2) = EV(:,:,2,1) >= EV(:,:,2,2);
            
                dist = compute.dist(V, V1, 4);
            
                if mod(iter_ct, 10) == 0
                    fprintf("\n\tIteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);
                end
            
                iter_ct = iter_ct + 1;
            
                V1 = V;            
            end

            fprintf("\n\tIteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);

        end

        function [mu1, Kagg] = getDist(g, amu, agrid, pil, pctDem)

            [nl, ~, ~, np] = size(g);
            nmu = max(size(amu));
            mu = zeros(nl, np, nmu); ixgrid = mu; wegrid = mu;

            for im = 1:nmu
                kval = amu(im);
                for il = 1:nl  
                    for ip = 1:np
                        % We linearly interpolate the policy function on bgrid to compute
                        % bfval on the distribution support.
                        
                        [ix, we] = compute.weight(agrid, kval);
                            
                        %split between rep and dem capital choices
                        kdval = pctDem*g(il,ix,1,ip)*we + ...
                                    pctDem*g(il,ix+1,1,ip)*(1.0 - we) + ...
                                    (1-pctDem)*g(il,ix,2,ip)*we + ...
                                    (1-pctDem)*g(il,ix+1,2,ip)*(1.0 - we);
                        
                        [ix, we] = compute.weight(amu, kdval);
                                
                        ixgrid(il, ip, im) = ix;
                        wegrid(il, ip, im) = we;
                    end
                end
            end
            
            distance = 20; iter_ct = 1;
            
            for i = 1:nmu
                for il = 1:nl
                    mu(il,1:np,i) = 1.0/(np*nmu*nl);
                end
            end

            while (distance > 1e-8)
                
                mu1 = zeros(size(mu));
                
                for im = 1:nmu
                    for il = 1:nl
                        for ip = 1:np
                            ix = ixgrid(il,ip,im); we = wegrid(il,ip,im);
                            muval = mu(il,ip,im);
                            
                            if muval > 0
                                for jl = 1:nl
                                    if (ix < nmu)
                                        mu1(jl,ip,ix) = mu1(jl,ip,ix) + pil(il,jl)*muval*we;
                                        mu1(jl,ip,ix+1) = mu1(jl,ip,ix+1) + pil(il,jl)*muval*(1.0 - we);
                                    else
                                        mu1(jl,ip,ix) = mu1(jl,ip,ix) + pil(il,jl)*muval*we;
                                    end
                                end
                            end
                        end                        
                    end
                end
                                
                distance = compute.dist(mu1, mu, 3);
                
                if (mod(iter_ct,25) == 0)
                    s = sprintf( '\n\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
                        iter_ct, distance, sum(sum(sum(mu1))));
                    disp (s);

%                     pause(1)
                    flatdistr=squeeze(sum(mu1,2));
                    mesh(flatdistr)
                end
                
                iter_ct = iter_ct + 1;
            
                mu = mu1;   
            end
            s = sprintf( '\n\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
                iter_ct, distance, sum(sum(sum(mu))));
            disp (s);

            distr2500 = sum(sum(mu,2),1);
            squeezed = squeeze(distr2500);
            Kagg = amu*squeezed;

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