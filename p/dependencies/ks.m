classdef ks
    methods(Static)

        function [Kprdata, Varray, Garray, EVarray] = getRegData(jt, terms, vTol, verbose)

            T = length(jt);
            pil = terms.pil;
            agrid = terms.agrid;
            lgrid = terms.lgrid;
            Kgrid = terms.Kgrid;
            na = length(agrid); nl = length(lgrid); nm = length(Kgrid);

            Kprdata = zeros(T,1);

%             gamma = terms.gamma;

            nmu = na*10;
            amu = linspace(agrid(1), agrid(na), nmu);
    
            fprintf("Solving HH problem...\n")
            terms.lambda = 0;
            [Varray, Garray, EVarray] = ks.parsolve(terms, vTol, verbose);

            Kprdata(1) = median(Kgrid);

            fprintf("Generating regression data...\n")
            for t = 2:1:T
                K_t = Kprdata(t-1);
                R = jt(t,1);
                d = jt(t,2);
                state_index = R*2-1 + d-1;


                % get the distribution

                [im, we] = compute.weight(Kgrid, K_t);

                G = Garray{state_index};
                g_int = we*G(im, :, :) + (1-we)*G(im + 1, :, :);
                g_int = squeeze(g_int);
                
                [~, Kagg] = HH.getDist(g_int, amu, agrid, pil, false);

                Kprdata(t) = Kagg;
                
                if mod(t,100) == 0
                    fprintf("\n\t t = %i", t)
                end
            end

        end


        function [Varray, Garray, EVarray] = parsolve(terms, vTol, verbose)

            Varray = cell(2,2);
            EVarray = cell(2,2);
            Garray = cell(2,2);
            RDmat = terms.RDmat;

            parfor i = 1:4
                j = ceil(i/2);
                k = mod(i,2) + 1;

                [V,G,EV] = ks.solve(terms, vTol, RDmat(i,:), j,k, verbose);

                Varray{i} = V;
                Garray{i} = G;
                EVarray{i} = EV;
            end

        end

        function [V, G, V0] = solve(terms, vTol, PIrd, Rind, dind, verbose)

            alpha = terms.alpha;
            phi = terms.phi;
            lgrid = terms.lgrid; nl = length(lgrid);
            agrid = terms.agrid; na = length(agrid);
            Kgrid = terms.Kgrid; nm = length(Kgrid);
            pil = terms.pil;
            g = terms.G;

            Kfore = terms.Kfore;
            eta = terms.etagrid(Rind);
            tau = terms.taugrid(Rind);
            terms.eta = eta;
            terms.tau = tau;
            delta = terms.dgrid(dind) + terms.delta;
            captax = terms.captax;

            state_index = Rind*2-1 + dind-1;

            V = zeros(nm, nl, na);
            G = V;
            V0 = V;
            TV = V; TG = G;

            % set up V so that it doesn't start empty
            % using believable r and w values: 4% interest, w = 1.3
            scale = .25;
            for im = 1:nm
                for ia = 1:na
                    kval = agrid(ia);
                    for il = 1:nl
                        yval = scale*(1+0.04*(1-captax(il)))*kval + ...
                            1.3*lgrid(il) - 0.04*phi;
                        ymin = max(1e-10, yval);
                        V(im, il, ia) = log(ymin);
                    end
                end
            end

            dist = 1e5;
            iter_ct = 1;
            
            %broadcasting, babyyyyyyyyyyy
            %forecasting K once, because the forecast only happens once
            
            Kpr = ks.forecast(Kfore, Kgrid', ...
                repmat(delta, nm, 1));

            %getting linear interpolation grid once too. 
            [ix, we] = ks.weight(Kgrid, Kpr);

            %then getting forecasted prices once, given future K.
            % (nR*ndelta) x nm grids
            rgrid = vaas.calcr(alpha, delta, Kpr, eta);
            wgrid = vaas.calcw(alpha, Kpr, eta);

            while dist > vTol
                
                EV = ks.getExpectation(V, ix, we, PIrd);
                for im = 1:nm
                    for ia = 1:na
                        for il = 1:nl
                            V0(im, il, ia) = pil(il,:)*EV(im, :,ia)';
                        end
                    end
                end

                for i = 1:nm
                    terms.r = rgrid(state_index, im);
                    terms.w = wgrid(state_index, im);
                    [TV(i,:,:), TG(i,:,:)]= spline.solve(terms,squeeze(V(i,:,:)), ...
                        squeeze(V0(i,:,:)), vTol);
                end 

                dist = compute.dist(V, TV, 3);
                kdist = compute.dist(G, TG, 3);
            
%                 if mod(iter_ct, 25) == 0
%                     fprintf("\n\tIteration %i: \n\t\t||TV - V|| = %4.6f" + ...
%                         "\n\t\t||TG - G|| = %4.6f", iter_ct, dist, kdist);
% %                     fprintf("\nInitial Values:");
% %                     fprintf("\nMin Kgrid: %2.4f, Max Kgrid: %2.4f", min(Kgrid), max(Kgrid));
% %                     fprintf("\nInitial Capital Forecasts: Kp = %2.4f, Kl = %2.4f", Kp, Kl);
% %                     fprintf("\nInitial Policy Function: Min Gp = %2.4f, Max Gp = %2.4f", ...
% %                         min(TGP(:)), max(TGP(:)));
% %                      fprintf("\nInitial Policy Function: Min Gl = %2.4f, Max Gl = %2.4f", ...
% %                         min(TGL(:)), max(TGL(:)));
% %                    fprintf("\nInitial Value Function: Min VP = %2.4f, Max VP = %2.4f", ...
% %                        min(TVP(:)), max(TVP(:)));
% %                    fprintf("\nInitial Value Function: Min VL = %2.4f, Max VL = %2.4f", ...
% %                        min(TVL(:)), max(TVL(:)));
%                 end
            
                iter_ct = iter_ct + 1;
            
                V = TV;
                G = TG;
            end

            if verbose
                fprintf("\n\tIteration %i (%i, %i): ||TV - V|| = %4.6f" + ...
                    "\t||TG - G|| = %4.6f\n", iter_ct, Rind, dind, dist, kdist);
            end

        end

        % get expectation
        % V is a three-dimensional object, because it's also defined on the
        % aggregate state grid. This function approximates expecations of V
        % given the aggregate state grid by linearly interpolating the
        % expecation over the aggregate states. The forecasting rule gives
        % our expected future state.
        function EV = getExpectation(V, ixmat, wemat, PIrd)

            [nm, ~, ~] = size(V);
            EV = zeros(size(V));
            Vproj = EV;

            for im = 1:nm
                for is = 1:length(PIrd)
                    ix = ixmat(is, im);
                    we = wemat(is, im);
                    Vproj = we.*V(ix, :, :) + (1-we).*V(ix+1, :, :);
                    EV(im,:,:) = EV(im,:,:) + PIrd(is)*Vproj;
                end
            end

        end

        function kpr = forecast(fore,k,d)
            % Ensure k and d are the same size for broadcasting
            if size(k,1) ~= size(d,1) || size(k,2) ~= size(d,2)
                error('K and d must have the same size.');
            end

            %k: nm x 1; d=nm x 1
            one_kd = [ones(size(k(:))) log(k(:)) log(d(:))]; % ==> one_kd = nm x 3
            %fore: (nR*nd) x 3
            kguess = fore*one_kd'; % ==> kguess = (nR*nd) x nm
            kpr =  exp(kguess);

        end


        function [ixmat, wemat] = weight(Kgrid, Kpr)

            Kpr = min(Kpr, ones(size(Kpr))*max(Kgrid));
            nm = length(Kgrid);

            % Kpr = (nR*nd) x nm
            % creating index matrix first

            ixmat = arrayfun(@(x) find(Kgrid <= x, 1, 'last'), Kpr);
            ixmat = min(ixmat, nm-1); % making sure it doesn't exceed

            Klower = Kgrid(ixmat); 
            Kupper = Kgrid(ixmat+1);

            wemat = (Kpr - Klower)./(Kupper - Klower);

        end

    end
end