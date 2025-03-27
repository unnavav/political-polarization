classdef predict
    methods(Static)

        %% simulate z 
        function inds = sim(N_T, Nz, rnseed, P_mat)
            sim = zeros(N_T, 1);
            isim = sim; %indices
        
            rng(rnseed);
            efsim = rand(N_T,1);
        
            %start at median
            isim(1) = round(median(1:Nz)); %should I round down? does median automatically?
            
            %build CDF
            cumP_mat = cumsum(P_mat, 2);
            
            %now, we simulate
            for t=1:N_T-1
                csumvec = cumP_mat(isim(t),1:Nz);
                condmet = (efsim(t+1) <= csumvec);
                isim(t+1) = find(condmet, 1, 'first');
            end
            
            inds = isim;
        end

        function [rt, wt, impliedK] = transition(r0, r1, lt, terms, lambda, dTol)
     
            verbose = false;

            alpha = terms.alpha;
            delta = terms.delta;

            agrid = terms.agrid;
            lgrid = terms.lgrid;
            pil = terms.pil;

            nl = length(lgrid); na = length(agrid);
            amu = linspace(min(agrid), max(agrid), na*10);

            T = length(lt);
            rt = linspace(r0, r1, T);
            wt = aiyagari.getW(rt, alpha, delta);
            
            impliedK = (rt + delta)/alpha;
            impliedK = impliedK.^(1/(alpha-1));
            impliedK = impliedK.*lt;
            Kguess = impliedK;

            % initalizing storage arrays
            Garray = cell(T,1);
            Varray = cell(T,1);
            Farray = cell(T,1);

            % solving final period value function
            terms.r = rt(T);
            terms.w = wt(T);
            [V,G,~,~] = HH.solve(nl, na, terms, 1e-8, verbose);
            Garray{T,1} = G;
            Varray{T,1} = V;
            [adistr, Kagg] = HH.getDist(G, amu, agrid, pil, true);
            Kguess(T) = Kagg;
            Farray{T,1} = adistr;

            % steady state in period 1
            terms.r = rt(1);
            terms.w = wt(1);
            [V,G,~,~] = HH.solve(nl, na, terms, 1e-8, verbose);
            Garray{1,1} = G;
            Varray{1,1} = V;
            [adistr, Kagg] = HH.getDist(G, amu, agrid, pil, true);
            Kguess(1) = Kagg;
            Farray{1,1} = adistr;

            DIST = 10;
            iter_ct = 1;

            while DIST > dTol
                 
                if DIST > 1
                    vTol = 1e-4;
                else
                    vTol = 1e-6;
                end
                fprintf("_Iteration %i_\n", iter_ct);
                for t = (T-1):-1:2
    
                    if mod(t,20) == 0
                        fprintf("T = %i... ", t)
                    end
                    
                    Vpr = Varray{t+1};
                    localterms = terms;
                    localterms.r = rt(t);
                    localterms.w = wt(t);
                    [V, G, ~] = HH.solve(nl, na, localterms, vTol, verbose);
    
                    Varray{t} = V;
                    Garray{t} = G;
                end
    
                for t = 2:T-1
                    
                    if mod(t,50) == 0
                        fprintf("T = %i...\n",t)
                    end
                    [Farray{t}, Kguess(t)] = HH.transitDistr(Garray{t}, ...
                        Farray{t-1}, amu, agrid, pil);
    
                end

                DIST = compute.dist(Kguess, impliedK,1);
                fprintf("||K - K'|| = %2.4f\n", DIST);

                rguess = alpha.*(Kguess./lt).^(alpha - 1) - delta;
                rguess(T) = rt(T);

                rt(2:T-1) = (1-lambda)*rt(2:T-1) + lambda*rguess(2:T-1);
                wt = aiyagari.getW(rt, alpha, delta);
                impliedK = (rt + delta)/alpha;
                impliedK = impliedK.^(1/(alpha-1));
                impliedK = impliedK.*lt;
                iter_ct = iter_ct + 1;
            end

            save transition_P_to_L_jan2025
        end

        %% perfect foresight

        function [zt, yt, ct, kt] = perfectForesight(zt, kst, k0, lambda, params,dTol)

            alpha = params(1);
            beta = params(2);
            delta = params(3);
            sigma = params(4);
            L = params(5);

            ist = delta*kstd;
            cst = kst^alpha - ist;

            T = size(zt, 2);
            zt = zt(1:T-1);
            kguess = repelem(kst, T-1);
            kguess(1) = k0*kst;
            
            kprguess = repelem(0, T-1);
            
            dist = max(abs(kguess - kprguess));
            
            iter_ct = 1;
            while dist > dTol
            
                yt = zt.*(kguess).^alpha*(L^(-alpha));
                rt = alpha.*zt.*(kguess.^(alpha-1))*(L^(1-alpha)) - delta;
            
                ct = repelem(cst,T-1);
                for i = flip(1:T-2)
                    if sigma == 1
                        ct(i) = ct(i+1)*((beta*(1+rt(i+1)))^(-1));
                    else
                        ct(i) = ct(i+1)*((beta*(1+rt(i+1)))^(-1/sigma));
                    end
                end
                ct(T-1) = cst;
            
                kt = repelem(kst, T);
                kt(1) = k0*kst;
                for i = 1:T-1
                    kt(i+1) = yt(i) - ct(i) + (1-delta)*kt(i);
                end
            
                kprguess = lambda*kguess + (1-lambda)*kt(1:T-1);
            
                dist = max(abs(kprguess - kguess)); 
            
                fprintf("Iteration %i: ||K' - K|| = %4.8f\n", iter_ct, dist);
            
                iter_ct = iter_ct + 1;
                kguess = kprguess;
            end
        end

    end
end