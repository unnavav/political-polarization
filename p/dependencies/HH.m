classdef HH
    methods(Static)


        function [Va, Vb, EVa, EVb, ga, gb] = solve(nl, na, np, terms, vTol, gTol)

            beta = terms.beta;
            sigma = terms.sigma;
            phi = terms.phi;
            r = terms.r;
            wage = terms.wage;
            lgrid = terms.lgrid;
            agrid = terms.agrid;
            pil = terms.pil;
            Y = terms.Y;
            p = terms.p;
            ubonus = terms.ubonus;

            Va = zeros(nl, na, np);
            EVa = Va;
            ga = Va;
            V1a = Va;
            Vb = Va;
            EVb = Vb;
            gb = Vb;
            V1b = Vb;

            for il = 1:nl
                parfor ia = 1:na
                    for ip = 1:np
                
                        max_cons = Y(il, ia, ip);
                         
                        if ip == 1
                            if sigma == 1
                                Va(il, ia, ip) = (log(max_cons)+ubonus)/(1-beta); 
                                Vb(il, ia, ip) = (log(max_cons))/(1-beta); 
                            else
                                Va(il, ia, ip) = ((max_cons^(1-sigma))/(1-sigma) + ...
                                    ubonus)/(1-beta);
                                Vb(il, ia, ip) = ((max_cons^(1-sigma))/(1-sigma))/(1-beta);
                            end
                        else
                            if sigma == 1
                                Vb(il, ia, ip) = (log(max_cons)+ubonus)/(1-beta); 
                                Va(il, ia, ip) = (log(max_cons))/(1-beta); 
                            else
                                Vb(il, ia, ip) = ((max_cons^(1-sigma))/(1-sigma) + ...
                                    ubonus)/(1-beta);
                                Va(il, ia, ip) = ((max_cons^(1-sigma))/(1-sigma))/(1-beta);
                            end
                        end
                    end
                end
            end

            dist = 10;
            iter_ct = 1;
            max_iter = 1000;
            while (dist > vTol && iter_ct <= max_iter)
                % get expected value function to start with
                for ip = 1:np
                    Vap1 = Va(:,:, 1); Vap2 = Va(:,:, 2);
                    Vbp1 = Vb(:,:, 1); Vbp2 = Vb(:,:, 2);
                    prob = p(ip);
                    parfor il = 1:nl
                        EVa(il, :, ip) = prob*pil(il,:)*Vap1 + ...
                            (1-prob)*pil(il,:)*Vap2;
                        EVb(il, :, ip) = prob*pil(il,:)*Vbp1 + ...
                            (1-prob)*pil(il,:)*Vbp2;
                    end
                end

                for il = 1:nl
                    for ip = 1:np
                    EVa_l = EVa(il,:,ip);
                    EVb_l = EVb(il,:,ip);
                        parfor ia = 1:na
                            y = Y(il, ia, ip);
                                    
                            if ip == 1
                                params = [y beta sigma ubonus];    
                                [res, aguess] = compute.gss(EVa_l, params, agrid, gTol);  
                                Va(il, ia, ip) = res;
                                ga(il, ia, ip) = aguess
                                params = [y beta sigma 0];    
                                [res, aguess] = compute.gss(EVb_l, params, agrid, gTol);  
                                Vb(il, ia, ip) = res;
                                gb(il, ia, ip) = aguess;
                            else
                                params = [y beta sigma 0];    
                                [res, aguess] = compute.gss(EVa_l, params, agrid, gTol);  
                                Va(il, ia, ip) = res;
                                ga(il, ia, ip) = aguess
                                params = [y beta sigma ubonus];    
                                [res, aguess] = compute.gss(EVb_l, params, agrid, gTol);  
                                Vb(il, ia, ip) = res;
                                gb(il, ia, ip) = aguess;
                            end
                        end
                    end
                end
                
                dist = max(compute.dist(Va, V1a, 3), ...
                    compute.dist(Vb, V1b, 3));
            
                if mod(iter_ct, 10) == 0
                    fprintf("\n\t\tIteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);
                end
            
                iter_ct = iter_ct + 1;
            
                V1a = Va;
                V1b = Vb;
            end
        
            fprintf("\n\t\tIteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);

        end

        function [mu1, Kagg] = getDist(ga, gb, amu, agrid, pil, pctA)

            [nl, ~, np] = size(ga);
            nmu = max(size(amu));
            mu = zeros(nl, nmu, np); 
            ixagrid = mu; weagrid = mu; ixbgrid = mu; webgrid = mu;
           
            for im = 1:nmu
                kval = amu(im);
                for il = 1:nl  
                    for ip = 1:np
                        % We linearly interpolate the policy function on bgrid to compute
                        % bfval on the distribution support.
                        
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
            
            for i = 1:nmu
                for il = 1:nl
                    for ip = 1:np
                        mu(il,i, ip) = 1.0/(nmu*nl*np);
                    end
                end
            end

            while (distance > 1e-8 && iter_ct < 3000)
                
                mu1 = zeros(size(mu));
                
                for im = 1:nmu
                    for il = 1:nl
                        for ip = 1:np

                            ixa = ixagrid(il,im, ip); wea = weagrid(il,im, ip); 
                            ixb = ixbgrid(il,im, ip); web = webgrid(il,im, ip); 
                            muval = mu(il, im, ip);
                            
                            if muval > 0
                                for jl = 1:nl
                                    if (ixa < nmu)
                                        mu1(jl,ixa, ip) = mu1(jl,ixa, ip) + ...
                                            pctA*pil(il,jl)*muval*we;
                                        mu1(jl,ixa+1, ip) = mu1(jl,ixa+1, ip) + ...
                                            pctA*pil(il,jl)*muval*(1.0 - we);
                                    else
                                        mu1(jl,ixa, ip) = mu1(jl,ixa,ip) + ...
                                            pctA*pil(il,jl)*muval*we;
                                    end
                                    
                                    if (ixb < nmu)
                                        mu1(jl,ixb, ip) = mu1(jl,ixb, ip) + ...
                                            (1-pctA)*pil(il,jl)*muval*we;
                                        mu1(jl,ixb+1, ip) = mu1(jl,ixb+1, ip) + ...
                                            (1-pctA)*pil(il,jl)*muval*(1.0 - we);
                                    else
                                        mu1(jl,ixb, ip) = mu1(jl,ixb,ip) + ...
                                            (1-pctA)*pil(il,jl)*muval*we;
                                    end
                                end
                            end

                        end
                    end
                end
                                
                distance = compute.dist(mu1, mu, 3);
                
%                 if (mod(iter_ct,200) == 0)
%                     s = sprintf( '\n\t\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
%                         iter_ct, distance, sum(sum(mu1)));
%                     disp(s);
%                 end
                
                iter_ct = iter_ct + 1;
            
                mu = mu1;   
            end
            s = sprintf( '\n\t\tIteration %3i: ||Tm - m|| = %8.6f\tsum = %6.4f ', ...
                iter_ct, distance, sum(sum(sum(mu))));
            disp(s);

            distr2500 = sum(sum(mu,3),1);
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