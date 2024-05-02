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