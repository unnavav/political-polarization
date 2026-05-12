classdef solve
    methods(Static)

        function [V, G, EV] = interpV(nl, na, EV, terms, vTol)

           dist = 10;
           iter_ct = 1;
           V = zeros(nl, na);
           G = zeros(nl, na);

           while dist > vTol && iter_ct < 500
               
                [TV, TG, EV] = compute.interpV(terms, V, EV, vTol);
        
                vdist = compute.dist(TV,V,2);
                gdist = compute.dist(TG,G,2);
                dist = max(vdist,gdist);
        
                if mod(iter_ct, 50) == 0
                    fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
                        "\n\t\t||TG - G|| = %4.6f", iter_ct, vdist, gdist);
                end
                iter_ct = iter_ct+1;
    
                V = TV;
                G = TG;
           end

        end

        function [kp, V, G] = getRegimeSS(nl, na, kval, K, EV, terms, vTol)

            kDist = 10;
            
            kp = [NA na];

            while kDist > dTol
    
                K = kval;
                fprintf("\n*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*-_-*")
                fprintf("\nA guess: %4.8f. Begin iteration for solution...\n", kval)
                fprintf("Guessing Probabilities:\n")
                
                pdist = 10;
                while pdist > dTol
                    EV = p*EVarray{1} + (1-p)*(EVarray{2});
                    % solve for VT:
                    terms.r = vaas.calcr(alpha, delta, K, epolicygrid(1));
                    terms.w = vaas.calcw(alpha, K, epolicygrid(1));
            
                    terms.tau = tpolicygrid(1);
                    % we have to get the value of lambda such that taxation 
                    % is redistributing everything. aka BB
                    tot_inc = (terms.w*lgrid).*ldist;
                    tot_inc = sum(tot_inc);
                    denom = (terms.w*lgrid).^(1-terms.tau).*ldist;
                    denom = sum(denom);
                    terms.lamval = tot_inc/denom;
            
                    [Varray{1}, Garray{1}] = solve.interpV(nl, na, EV, terms, dTol);
            
                    terms.r = vaas.calcr(alpha, delta, K, epolicygrid(2));
                    terms.w = vaas.calcw(alpha, K, epolicygrid(2));
            
                    terms.tau = tpolicygrid(2);
                    % we have to get the value of lambda such that taxation 
                    % is redistributing everything. aka BB
                    tot_inc = (terms.w*lgrid).*ldist;
                    tot_inc = sum(tot_inc);
                    denom = (terms.w*lgrid).^(1-terms.tau).*ldist;
                    denom = sum(denom);
                    terms.lamval = tot_inc/denom;
            
                    [Varray{2}, Garray{2}] = solve.interpV(nl, na, EV, terms, dTol);
            
                    EVarray{1} = predict.getExpectation(Varray{1}, pil);
                    EVarray{2} = predict.getExpectation(Varray{2}, pil);
            
                    [Warray{1}, Karray{1}] = HH.getDist(Garray{1}, amu, agrid, ...
                        pil, false);  
            
                    acond = compute.condense(Warray{1}, amu, agrid);
                    VOTES = EVarray{1,1} > EVarray{1,2};
                    VOTES = acond .* VOTES;
                    pnew = sum(sum(VOTES));
                    
                    pdist = abs(pnew-p);
                    fprintf("\tProbability diff: |%0.2f - %0.2f| = %0.2f\n", pnew, ...
                        p, pdist);
                    p = (p+pnew)/2;
                end
            
                kp = [kp; kval p];
            
                [Warray{2}, Karray{2}] = HH.getDist(Garray{2}, amu, agrid, ...
                      pil, false);            
                kdist = Karray{1} - kval;
            
                if kdist > 0
                    fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too low.\n", abs(kdist))
                    kl = (1-adj)*kval+adj*kl;
                else
                    fprintf("\n||Kguess - Kagg|| = %4.5f. \tAggregate capital is too high.\n", abs(kdist))
                    kh = (1-adj)*kval+(adj)*kh;
                end
            
                kDist = abs(kdist);  
                kval = .5*(kl + kh);
            end

        end
    end
end