classdef predict
    methods(Static)

        %% simulate z 
        function zinds = simz(N_T, Nz, rnseed, P_mat)
            zsim = zeros(1,N_T);
            izsim = zsim; %indices
        
            rng(rnseed);
            efsim = rand(1,N_T);
        
            %start at median
            izsim(1) = median(1:Nz); %should I round down? does median automatically?
            
            %build CDF
            cumP_mat = cumsum(P_mat, 2);
            
            %now, we simulate
            for t=1:N_T-1
                csumvec = cumP_mat(izsim(t),1:Nz);
                condmet = (efsim(t+1) <= csumvec);
                izsim(t+1) = find(condmet, 1, 'first');
            end
            
            zinds = izsim;
        end

        function [rt] = transition(r0, r1, lt, terms, lambda, dTol)
     
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
            wt = aiyagari.getW(rtguess, alpha, delta);
            
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
            [G,V,~,~] = HH.solve(nl, na, terms, dTol, verbose);
            Garray{T,1} = G;
            Varray{T,1} = V;
            [adistr, Kagg] = HH.getDist(G, amu, agrid, pil, true);
            Kguess(T) = Kagg;
            Farray{T,1} = adistr;

            % steady state in period 1
            terms.r = rt(1);
            terms.w = wt(1);
            [G,V,~,~] = HH.solve(nl, na, terms, dTol, verbose);
            Garray{1,1} = G;
            Varray{1,1} = V;
            [adistr, Kagg] = HH.getDist(G, amu, agrid, pil, true);
            Kguess(1) = Kagg;
            Farray{1,1} = adistr;


            DIST = 10;
            iter_ct = 1;

            while DIST > dTol
                 
                fprintf("_Iteration %i_\n", iter_ct);
                for t = (T-1):-1:2
    
%                     if mod(t,50) == 0
%                         disp(t)
%                     end
                    Vpr = Varray{t+1};
                    localterms = terms;
                    localterms.r = rt(t);
                    localterms.w = wt(t);
                    [V, G,~] = HH.solve(nl, na, localterms, dTol/10, verbose);
    
                    Varray{t} = V;
                    Garray{t} = G;
                    mesh(G);
                end
    
                for t = 2:T-1
                    
%                     if mod(t,50) == 0
%                         disp(t)
%                     end
                    [Farray{t}, Kguess(t)] = HH.transitDistr(Garray{t}, ...
                        Farray{t-1}, amu, agrid, pil);
    
                end

                DIST = compute.dist(Kguess, impliedK,1);
                fprintf("||K - K'|| = %2.4f\n", DIST);

                rguess = alpha.*(Kguess./lt).^(alpha - 1) - delta;

                rt = (1-lambda)*rt + lambda*rguess;
                impliedK = (rt + delta)/alpha;
                impliedK = impliedK.^(1/(alpha-1));
                impliedK = impliedK.*lt;
                iter_ct = iter_ct + 1;
            end

        end

        %% perfect foresight

        function [zt, yt, ct, kt] = perfectForesight(zt, kst, k0, lambda, params,dTol)

            alpha = params(1);
            beta = params(2);
            delta = params(3);
            sigma = params(4);
            L = params(5);

            ist = delta*kst;
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