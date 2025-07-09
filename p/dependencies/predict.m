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

            terms.phi = 0;

            % solving final period value function
            terms.r = rt(T);
            terms.w = wt(T);
            [V,G] = predict.solve(nl, na, terms, agrid, lgrid);
            Garray{T,1} = G;
            Varray{T,1} = V;
            [adistr, Kagg] = HH.getDist(G, amu, agrid, pil, true);
            Kguess(T) = Kagg;
            Farray{T,1} = adistr;

            % steady state in period 1
            terms.r = rt(1);
            terms.w = wt(1);
            [V,G] = predict.solve(nl, na, terms, agrid, lgrid);
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
    
%                     if mod(t,20) == 0
                        fprintf("T = %i... ", t)
%                     end
                    
                    Vpr = Varray{t+1};

                    localterms = terms;
                    localterms.r = rt(t);
                    localterms.w = wt(t);
               
                    
                    [V, G] = predict.solve(nl,na,localterms,agrid,lgrid);

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

        %% perfect foresight for regimes

        function [Rguess, pt, EVarray, Garray, Warray] = perfectForesight(kt, p0, terms, dTol)

            T = length(kt);
            etagrid = terms.etagrid;
            taugrid = terms.taugrid;
            agrid = terms.agrid;
            amu = terms.amu;
            alpha = terms.alpha;
            delta = terms.delta;
            pil = terms.pil;

            nl = length(terms.lgrid); na = length(terms.agrid);

            Varray = cell(T, 2);
            Garray = cell(T, 2);
            EVarray = cell(T,2);
            Warray = cell(T,1);

            %dynamic tracking
            ptguess = ones(T,1)*p0;
            Rguess = ptguess<.5+1;
            pt = ptguess;
            Rt = Rguess;
            EVgap = zeros(T,1);

            fprintf("Solving end steadystate vals.\n")
            EVguess = zeros(length(terms.lgrid), length(terms.agrid));
            K = kt(T);
            % solve for VT:
            terms.r = vaas.calcr(alpha, delta, K, etagrid(1));
            terms.w = vaas.calcw(alpha, K, etagrid(1));
            terms.tau = taugrid(1);
            [Varray{T,1}, Garray{T,1}] = solve.interpV(nl, na, EVguess, terms, dTol/100);

            terms.r = vaas.calcr(alpha, delta, K, etagrid(2));
            terms.w = vaas.calcw(alpha, K, etagrid(2));
            terms.tau = taugrid(2);
            [Varray{T,2}, Garray{T,2}] = solve.interpV(nl, na, EVguess, terms,  dTol/100);

            EVarray{T,1} = predict.getExpectation(Varray{T,1}, pil);
            EVarray{T,2} = predict.getExpectation(Varray{T,2}, pil);
            [Warray{T}, ~] = HH.getDist(Garray{T,(p0>.5)+1}, amu, agrid, pil, false);
            acond = compute.condense(Warray{T}, amu, agrid);
            VOTES = EVarray{T,1} > EVarray{T,2};
            VOTES = acond .* VOTES;
            ptguess(T) = sum(sum(VOTES));
            Rguess(T) = (ptguess(T) > 0.5) + 1;

            Rguess(T) = (p0>.5)+1;

            fprintf("\nSolving beginning steadystate vals.\n")
            % period 1 value function
            K = kt(1);
            % solve for VT:
            terms.r = vaas.calcr(alpha, delta, K, etagrid(1));
            terms.w = vaas.calcw(alpha, K, etagrid(1));
            terms.tau = taugrid(1);
            [Varray{1,1}, Garray{1,1}] = solve.interpV(nl, na, EVguess, terms, dTol);

            terms.r = vaas.calcr(alpha, delta, K, etagrid(2));
            terms.w = vaas.calcw(alpha, K, etagrid(2));
            terms.tau = taugrid(2);
            [Varray{1,2}, Garray{1,2}] = solve.interpV(nl, na, EVguess, terms, dTol);

            EVarray{1,1} = predict.getExpectation(Varray{1,1}, pil);
            EVarray{1,2} = predict.getExpectation(Varray{1,2}, pil);

            [Warray{1}, ~] = HH.getDist(Garray{1,Rt(1)}, amu, agrid, pil, false);
            acond = compute.condense(Warray{1}, amu, agrid);
            VOTES = EVarray{1,1} > EVarray{1,2};
            VOTES = acond .* VOTES;
            ptguess(1) = sum(sum(VOTES));
            Rguess(1) = (ptguess(1) > 0.5) + 1;

            %to do: 1. make sure I get starting and ending regimes right
            % 2. make sure I'm getting regime convergence right
            % 3. save the distributions so I can solve the distributions 
            % forward.

            pdist = 10;
            iter = 1;
            while pdist > dTol

                fprintf("\nIteration %i \n", iter);
                fprintf("\tBackwards pass...\n")
                % backwards pass
                for t = T-1:-1:1
                    if mod(t,50) == 0
                        fprintf("t = %i ...", t);
                    end
                    K = kt(t);
                    pguess = pt(t+1);
                    pnew = 0; inner_iter = 1; dist = 10;
                    EV = pguess*EVarray{t+1,1} + (1-pguess)*(EVarray{t+1,2});
                    
                    while dist > 0.001
                        R = (pguess < .5) + 1;        
                        terms.r = vaas.calcr(alpha, delta, K, etagrid(1));
                        terms.w = vaas.calcw(alpha, K, etagrid(1));
                        terms.tau = taugrid(1);
                        [Varray{t,1}, Garray{t,1}] = compute.backsolve(terms, EV, dTol);
    
                        terms.r = vaas.calcr(alpha, delta, K, etagrid(2));
                        terms.w = vaas.calcw(alpha, K, etagrid(2));
                        terms.tau = taugrid(2);
                        [Varray{t,2}, Garray{t,2}] = compute.backsolve(terms, EV, dTol);
    
                        EVarray{t,1} = predict.getExpectation(Varray{t,1}, pil);
                        EVarray{t,2} = predict.getExpectation(Varray{t,2}, pil);
                        EVgap(t) = sum(sum(EVarray{t,1}-EVarray{t,2}));
                        EV = pguess*EVarray{t+1,1} + (1-pguess)*(EVarray{t+1,2});

                        [W, K] = HH.getDist(Garray{1,R}, amu, agrid, pil, false);   
                        acond = compute.condense(W, amu, agrid);
                        VOTES = EVarray{t,1} > EVarray{t,2};
                        VOTES = acond .* VOTES;
                        pnew = sum(sum(VOTES));
                        if inner_iter > 100
                            warning("pt(t) not converging at t = %d", t);
                            break
                        elseif mod(inner_iter,50) == 0
                            fprintf("\t\tInner iter: %i \n", inner_iter)
                        end
                        inner_iter=inner_iter+1;
                        dist = abs(pguess-pnew);
                        pguess = pnew;
                    end
                    pt(t) = pnew;
                end
                Rt = (pt<.5)+1;
                
                [Warray{1}, ~] = HH.getDist(Garray{1,Rt(1)}, amu, agrid, pil, false);   
                acond = compute.condense(Warray{1}, amu, agrid);
                VOTES = EVarray{1,1} > EVarray{1,2};
                VOTES = acond .* VOTES;
                ptguess(1) = sum(sum(VOTES));
                Rguess(1) = (ptguess(1) > 0.5) + 1;

                fprintf("\tForwards pass...\n")
                % forwards pass: updating capital path
                for t = 2:1:T-1
                    if mod(t,50) == 0
                        fprintf("t = %i ...", t);
                    end  
                    % get today's distribution
                    G = Garray{t};
                    [Warray{t}, Kagg] = HH.transitDistr(G, Warray{t-1}, ...
                        amu, agrid, pil);
    
                    acond = compute.condense(Warray{t,1}, amu, agrid);
                    VOTES = EVarray{t,1} > EVarray{t,2};
                    VOTES = acond.*VOTES;
                    ptguess(t) = sum(sum(VOTES));
    
                    Rguess(t) = (ptguess(t)>.5) + 1;
                end
                
                pdist = sum(abs(ptguess - pt));  
                pt = .5*ptguess + .5*pt;
                fprintf("\n Pdist = %0.4f \n", pdist);
                iter = iter+1;
            end

        end

        % get expectation 
        function EV = getExpectation(V, pil)

            [nl, na] = size(V);
            EV = zeros(size(V));
            for ia = 1:na
                for il = 1:nl
                    EV(il, ia) = pil(il,:)*V(:,ia);
                end
            end
        end
    end
end