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

        function [rt, Yt, Ct, Kt] = perfectForesight(rstart, rend, L, lambda, terms, dTol)

            alpha = terms.alpha;
            beta = terms.beta;
            delta = terms.delta;
            sigma = terms.sigma;

            T = size(L, 2);
            rt = zeros(1,T);
            Yt = rt; Ct = rt; It = rt; Kt = rt;

            r0guess = linspace(rstart, rend, T-1);

            Kt = ((r0guess + delta)/alpha).^(1/(alpha-1)).*L;
            Ks = ((rstart + delta)/alpha).^(1/(alpha-1)).*L;
            
            dist = max(abs(r1guess - r0guess));
            
            iter_ct = 1;

            while dist > dTol

                % given prices, calculate household problem
                % so first get this period's wages given r_t, L_t

                Yt(T) = Kt(T)^alpha * L(T)^(1-alpha);
                Ct(T) = Yt(T) - delta*Kt(T);

                for i = flip(1:T-2)
                    Yt(i) = Kt(i)^alpha * L(i)^(1-alpha);
                    %keeping homogenous capital for now
                    Ct(i) = max((beta * (1 / Ct(i+1) * (1 + r0guess(i+1) * (1 - 0.15))))^(1 / sigma), 1e-8);                    
                end

                % now get implied K
                for i = 1:T-1
                    Ks(i+1) = Yt(i) - Ct(i) + (1-delta)*Ks(i);
                end

                kprguess = lambda*Kt(1:T-1) + (1-lambda)*Ks(1:T-1);
            
                dist = max(abs(alpha*(Ks(T)/L(T))^(alpha - 1) - delta - rend)); 
            
                fprintf("Iteration %i: ||K' - K|| = %4.8f\n", iter_ct, dist);
            
                iter_ct = iter_ct + 1;
                Kt = [kprguess kend];
                r0guess = alpha*(Kt./L).^(alpha - 1) - delta;
            end

            rt = alpha*(Kt./L).^(alpha - 1) - delta;
        end

    end
end