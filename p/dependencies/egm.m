classdef egm
    methods(Static)

        function [TV, G, C, V0] = solve(pol_terms, V0, V)
            
            [nl, na] = size(V);

            w = pol_terms.w;
            r = pol_terms.r;

            captax = pol_terms.captax;
            lamval = pol_terms.lamval;
            tau = pol_terms.tau;
            agrid = pol_terms.agrid;
            lgrid = pol_terms.lgrid;
            beta = pol_terms.beta;
            sigma = pol_terms.sigma;

            phi = 1;

            endoK = zeros(size(V));
            TV = endoK;
            D = egm.numdev(V0, agrid);
            for ia = 1:na
                kpr = agrid(ia);
                for il = 1:nl
                    l = lgrid(il);

                    c = (beta*D(il, ia))^(-1/sigma);
                    y = gov.tax(w*l, lamval, tau);
                    
                    numerator = c + kpr - y;
                    denom = (1+r*(1-captax(il))); % TODO: check phi in right place

                    impliedK = numerator/denom;

                    endoK(il, ia) = impliedK;
                end
            end
        
            G = zeros(size(V));
            %linterpolate decision rule
            for ia = 1:na
                k = agrid(ia);
                for il = 1:nl
                    lkvals = endoK(il, :);
                    
                    if k <= lkvals(1)
                        G(il, ia) = agrid(1);
                    elseif k > lkvals(end)
                        G(il,ia) = agrid(na);
                    else 
                        [ix, we] = compute.weight(lkvals, k);
                        kpr = we*agrid(ix) + (1-we)*agrid(ix + 1);
                        G(il, ia) = max(agrid(1), kpr);
                    end
                end
            end
            G = max(agrid(1), G);
            G = min(agrid(na),G);
        
            C = endoK;
            for ia = 1:na
                for il = 1:nl
                    l = lgrid(il);

                    c = (1 + r*(1 - captax(il))) * agrid(ia) ...
                        + gov.tax(w*l, lamval, tau) - G(il, ia);

                    C(il, ia) = max(1e-6, c);

                    [ix, we] = compute.weight(agrid, G(il, ia));
                    ev = we*V0(il, ix) + (1-we)*V0(il, ix+1);
                    TV(il, ia) = HH.u(C(il,ia), sigma) + beta*ev;
                end
            end
        end

        %% getting derivative of expected value function
        % inputs: EV grid
        % outputs: DEV grid
        function DEV = numdev(V0, agrid)
            [~, na] = size(V0);
            
            DEV = zeros(size(V0));

            for ia = 1:na
                evia = V0(:, ia);
                a = agrid(ia);

                if (ia > 1 && ia < na)
                    evia_0 = V0(:, ia-1);
                    evia_1 = V0(:, ia+1);
                    a1 = agrid(ia+1);
                    a0 = agrid(ia - 1);
                    
                    dr = (evia_1 - evia)/(a1 -a);
                    dl = (evia - evia_0)/(a -a0);

                    % scaling by spacing rather than taking a simple avg
                    DEV(:,ia) = ((a - a0)/(a1 - a0)).*dr + ...
                        ((a1 - a)/(a1 - a0)).*dl;

                elseif ia == 1

                    evia_1 = V0(:, ia+1);
                    a1 = agrid(ia+1);
                    dr = (evia_1 - evia)/(a1 -a);

                    DEV(:, ia) = dr;

                else

                    evia_0 = V0(:, ia-1);
                    a0 = agrid(ia-1);
                    dl = (evia - evia_0)/(a - a0);

                    DEV(:, ia) = dl;

                end

            end

            DEV = max(min(DEV, 1e12), 1e-12);
                    
        end

        function Dev = aubdev(ev,agrid)
            [m, anum] = size(ev);

            Dev = zeros(m, anum);

            for ia = 1:1:anum
        
                % numerical derivative as secant
                if (ia == anum)
                    % left hand limit if afval is upper bound for agrid
                    Dev(1:m, ia) = Devr(1:m,1); % (ev(1:m,ia) - ev(1:m,ia-1))./(agrid(ia) - agrid(ia-1));
                elseif(ia == 1)
                    % right hand limit if afval is less than upper bound
                    Dev(1:m,ia) = (ev(1:m,ia+1) - ev(1:m,ia))./(agrid(ia+1) - agrid(ia));
                    Devl(1:m,1) = Dev(1:m,ia);
                else
                    %Devl(1:m,ia) = (ev(1:m,ia) - ev(1:m,ia-1))./(agrid(ia) - agrid(ia-1));
                    Devr(1:m,1) = (ev(1:m,ia+1) - ev(1:m,ia))./(agrid(ia+1) - agrid(ia));
                    Dev(1:m,ia) = (Devl + Devr)/2.0;
                    Devl = Devr;
                end
        
            end
        end

       function d = solveD(EV0, ik, kchgrid)
            nk = length(kchgrid);
            k = kchgrid(ik);

            if ik == nk
                d = (EV0(ik) - EV0(ik-1))/(k - kchgrid(ik-1));
            else
                d = (EV0(ik+1) - EV0(ik))/(kchgrid(ik+1) - k);
            end
        end

    end
end
