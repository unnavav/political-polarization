classdef ks
    methods(Static)

        function [Kdata, adist] = getRegData(jt, terms, foreguess, nK, dTol, vTol, verbose)

            T = length(jt);
            agrid = terms.agrid;
            lgrid = terms.lgrid;
            na = length(agrid); nl = length(lgrid);
            Kgrid = linspace(agrid(1), agrid(na), nK);

            Kdata = zeros(1,T);

            for t = 1:1:T

                while DIST > dTol

                    % nm is the number of states (?)
                    K = [repelem(1, nm); log(Kguess)];
                    Kfore = exp(foreguess*K);
        
                    [V, G, C, V0] = ks.solve(nl, na, terms, vTol, verbose);
                end
            end

        end

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

            EV = ks.getExpectation(V,Kgrid,)
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


    end
end