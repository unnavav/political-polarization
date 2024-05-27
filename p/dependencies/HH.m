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
            g = terms.G;
            p = terms.p;
            ubonus = terms.ubonus;
            captax = terms.captax;

            Va = zeros(nl, na, np);
            EVa = Va;
            ga = Va;
            V1a = Va;
            Vb = Va;
            EVb = Vb;
            gb = Vb;
            V1b = Vb;

            for il = 1:nl
                for ia = 1:na
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
                    for il = 1:nl
                        EVa(il, :, ip) = prob*pil(il,:)*Vap1 + ...
                            (1-prob)*pil(il,:)*Vap2;
                        EVb(il, :, ip) = prob*pil(il,:)*Vbp1 + ...
                            (1-prob)*pil(il,:)*Vbp2;
                    end
                end

                % commence EGM 😎😎

                %step 1: get derivatives

                DEVa = egm.numdev(EVa, agrid);
                DEVb = egm.numdev(EVb, agrid);

                %step 2: get endogenous grid
                Ea = zeros(size(DEVa));
                Eb = zeros(size(DEVa));
                for il = 1:nl
                    for ip = 1:np
                        prob = p(ip);
                        for ia = 1:na
                            deva = prob*DEVa(il, ia, 1) + (1-prob)*DEVa(il, ia, 2);
                            devb = prob*DEVb(il, ia, 1) + (1-prob)*DEVb(il, ia, 2);
                            apr = agrid(ia);
                            if sigma == 1
                                ca = exp(beta*deva);
                                cb = exp(beta*devb);
                            else
                                cb = (beta*deva)^(-sigma);
                                ca = (beta*devb)^(-sigma);
                            end

                            Ea(il, ia, ip) = (ca + apr - Y(il, ia, ip) - ...
                                g(ip))/(1+r*(1-captax)); 
                            Eb(il, ia, ip) = (cb + apr - Y(il, ia, ip) - ...
                                g(ip))/(1+r*(1-captax)); 
                        end
                    end
                end

                %step 3: back out normal vals
                for il = 1:nl
                    for ip = 1:np
                        lb_a = min(Ea(il, :, ip));
                        lb_b = min(Ea(il, :, ip));

                        for ia = 1:na
                            ahat = agrid(ia);
            
                            if ahat < lb_a
                                ga(il, ia, ip) = 0;
                            else
                                [ix, we] = compute.weight(Ea(il, :, ip),ahat);
                                ga(il, ia, ip) = we*Ea(il, ix, ip) + (1-we)*Ea(il, ix+1, ip);
                            end

                            if ahat < lb_b
                                gb(il, ia, ip) = 0;
                            else
                                [ix, we] = compute.weight(Eb(il, :, ip), ahat);
                                gb(il, ia, ip) = we*Eb(il, ix, ip) + (1-we)*Eb(il, ix+1, ip);
                            end
            
                            c_a = (1+r)*ahat + Y(il, ia, ip) + g(ip) - r*phi - ga(il, ia, ip);
                            c_b = (1+r)*ahat + Y(il, ia, ip) + g(ip) - r*phi - gb(il, ia, ip);
            
                            if (c_a < 0 || c_b < 0)
                                disp([il ia ip])
                                return
                            end

                            if ip == 1
                                bonusa = ubonus;
                                bonusb = 0;
                            else
                                bonusa = 0;
                                bonusb = ubonus;
                            end

                            if sigma == 1
                                Va(il, ia, ip) = log(c_a) + bonusa + beta*EVa(il, ia, ip);
                                Vb(il, ia, ip) = log(c_b) + bonusb + beta*EVb(il, ia, ip);
                            else
                                Va(il, ia, ip) = (c_a^(1-sigma))/(1-sigma) + ...
                                    bonusa + beta*EVa(il, ia, ip);
                                Vb(il, ia, ip) = (c_b^(1-sigma))/(1-sigma) + ...
                                    bonusb + beta*EVb(il, ia, ip);
                            end
                        end
                    end
                end
                 

%                 for il = 1:nl
%                     for ip = 1:np
%                     EVa_l = EVa(il,:,ip);
%                     EVb_l = EVb(il,:,ip);
%                         parfor ia = 1:na
%                             y = Y(il, ia, ip);
%                                     
%                             if ip == 1
%                                 params = [y beta sigma ubonus];    
%                                 [res, aguess] = compute.gss(EVa_l, params, agrid, gTol);  
%                                 Va(il, ia, ip) = res;
%                                 ga(il, ia, ip) = aguess
%                                 params = [y beta sigma 0];    
%                                 [res, aguess] = compute.gss(EVb_l, params, agrid, gTol);  
%                                 Vb(il, ia, ip) = res;
%                                 gb(il, ia, ip) = aguess;
%                             else
%                                 params = [y beta sigma 0];    
%                                 [res, aguess] = compute.gss(EVa_l, params, agrid, gTol);  
%                                 Va(il, ia, ip) = res;
%                                 ga(il, ia, ip) = aguess
%                                 params = [y beta sigma ubonus];    
%                                 [res, aguess] = compute.gss(EVb_l, params, agrid, gTol);  
%                                 Vb(il, ia, ip) = res;
%                                 gb(il, ia, ip) = aguess;
%                             end
%                         end
%                     end
%                 end
                
                dist = max(compute.dist(Va, V1a, 3), ...
                    compute.dist(Vb, V1b, 3));
            
                if mod(iter_ct, 100) == 0
                    fprintf("\n\t\tIteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);
                end
            
                iter_ct = iter_ct + 1;
            
                V1a = Va;
                V1b = Vb;
            end
        
            fprintf("\n\t\tIteration %i: ||TV - V|| = %4.7f\n", iter_ct, dist);

        end

        function [mu1A, KaggA, mu1B, KaggB] = getDist(ga, gb, amu, agrid, pil, pctA)

            [nl, ~, np] = size(ga);
            nmu = max(size(amu));
            muA = zeros(nl, nmu);
            muB = zeros(nl, nmu);
            ixagrid = zeros(size(ga)); weagrid = ixagrid; 
            ixbgrid = ixagrid; webgrid = ixagrid;
           
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
            
            for im = 1:nmu
                for il = 1:nl
                    for ip = 1:np
                        muA(il,im) = 1.0/(nmu*nl);
                        muB(il,im) = 1.0/(nmu*nl);
                    end
                end
            end

            while (distance > 1e-6 && iter_ct < 3000)
                
                mu1A = zeros(size(muA));
                mu1B = zeros(size(muB));
                
                for im = 1:nmu
                    for il = 1:nl
                        ip = 1;

                        ixa = ixagrid(il,im, ip); wea = weagrid(il,im, ip); 
                        ixb = ixbgrid(il,im, ip); web = webgrid(il,im, ip); 
                        muvalA = muA(il, im); 
                        
                        %first do A distr
                        if muvalA > 0
                            for jl = 1:nl
                                %a households' movements
                                if (ixa < nmu)
                                    mu1A(jl,ixa) = mu1A(jl,ixa) + ...
                                        pctA*pil(il,jl)*muvalA*wea;
                                    mu1A(jl,ixa+1) = mu1A(jl,ixa+1) + ...
                                        pctA*pil(il,jl)*muvalA*(1.0 - wea);

                                else
                                    mu1A(jl,ixa) = mu1(jl,ixa) + ...
                                        pctA*pil(il,jl)*muvalA;
                                    mu1A(jl,ixa) = mu1(jl,ixa) + ...
                                        pctA*pil(il,jl)*muvalA;
                                end
                                
                                %b households' movements
                                if (ixb < nmu)
                                    mu1A(jl,ixb) = mu1A(jl,ixb) + ...
                                        (1-pctA)*pil(il,jl)*muvalA*web;
                                    mu1A(jl,ixb+1) = mu1A(jl,ixb+1) + ...
                                        (1-pctA)*pil(il,jl)*muvalA*(1.0 - web);
                                else
                                    mu1A(jl,ixb) = mu1A(jl,ixb) + ...
                                        (1-pctA)*pil(il,jl)*muvalA;
                                end
                            end
                        end

                        % now B distr
                        ip = 2;
                        ixa = ixagrid(il,im, ip); wea = weagrid(il,im, ip); 
                        ixb = ixbgrid(il,im, ip); web = webgrid(il,im, ip);
                        muvalB = muB(il, im);

                        if muvalB > 0
                            for jl = 1:nl
                                % a HH movements for B in power
                                if (ixa < nmu)
                                    mu1B(jl,ixa) = mu1B(jl,ixa) + ...
                                        pctA*pil(il,jl)*muvalB*wea;
                                    mu1B(jl,ixa+1) = mu1B(jl,ixa+1) + ...
                                        pctA*pil(il,jl)*muvalB*(1.0 - wea);
                                else
                                    mu1B(jl,ixa) = mu1B(jl,ixa) + ...
                                        pctA*pil(il,jl)*muvalB;
                                    mu1B(jl,ixa) = mu1B(jl,ixa) + ...
                                        pctA*pil(il,jl)*muvalB;
                                end
                                
                                if (ixb < nmu)
                                    mu1B(jl,ixb) = mu1B(jl,ixb) + ...
                                        (1-pctA)*pil(il,jl)*muvalB*web;
                                    mu1B(jl,ixb+1) = mu1B(jl,ixb+1) + ...
                                        pctA*pil(il,jl)*muvalB*(1.0 - web);
                                else
                                    mu1B(jl,ixb) = mu1B(jl,ixb) + ...
                                        (1-pctA)*pil(il,jl)*muvalB*we;
                                end
                            end
                        end

                    end
                end
                                
                distanceA = compute.dist(mu1A, muA, 3);
                distanceB = compute.dist(mu1B, muB, 3);

                if (mod(iter_ct,200) == 0)
                    s = sprintf( '\n\t\tIteration %3i: ||Tma - ma|| = %8.6f\tsum = %6.4f ', ...
                        iter_ct, distanceA, sum(sum(mu1A)));
                    disp(s);
                    s = sprintf( '\n\t\tIteration %3i: ||Tmb - mb|| = %8.6f\tsum = %6.4f ', ...
                        iter_ct, distanceB, sum(sum(mu1B)));
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