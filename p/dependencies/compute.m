classdef compute
    methods(Static)
        function [res, aval, v] = gss(Cvals, y, beta, sigma, searchgrid, prec)            
            ni = max(size(searchgrid));
            r = (3-sqrt(5))/2;

            a = searchgrid(1);
            b = min(searchgrid(ni), y); % Ensure c > 0;
            c = (1-r)*a + r*b;
            d = r*a + (1-r)*b;

            vc = compute.linterpolate(Cvals, searchgrid, c);
            fc = -HH.util(c, y, beta, sigma, vc);

            vd = compute.linterpolate(Cvals, searchgrid, d);
            fd = -HH.util(d, y, beta, sigma, vd);
            iter_ct = 1;

%             disp(fc);
%             disp(fd);
            v = [c fc; d fd];

            while abs(a - b) > prec
                if fc > fd
                    a = c;
                    c = d;
                    d = r*a + (1-r)*b;

                    fc = fd;
                    vd = compute.linterpolate(Cvals, searchgrid, d);
                    fd = -HH.util(d, y, beta, sigma, vd);
                    v = [v; d fd];
                else
                    b = d;
                    d = c;
                    c = (1-r)*a + r*b;

                    fd = fc;
                    vc = compute.linterpolate(Cvals, searchgrid, c);
                    fc = -HH.util(c, y, beta, sigma, vc);
                    v = [v; c fc];
                end

%                 fprintf("Results of iteration %i: [a c d b] = [%4.4f %4.4f %4.4f %4.4f]\n", ...
%                     iter_ct, a, c, d, b);
% 
%                 disp(fc);

                iter_ct = iter_ct+1;
            end

%             disp(fc);

            aval = c;
            vc = compute.linterpolate(Cvals, searchgrid, c);
            res = HH.util(c, y, beta, sigma, vc);
        end

        function v = linterpolate(Vvec, searchgrid, vi)
            ni = max(size(searchgrid));
            if vi <= searchgrid(1)
%                 il = 1;
                v = Vvec(1);
            elseif vi >= searchgrid(ni)
%                 il = nk-1;
                v = Vvec(ni);
            else
                il = sum(searchgrid <= vi);
                wl = (searchgrid(il + 1) - vi)/(searchgrid(il+1)-searchgrid(il));
                v = wl*Vvec(il) + (1-wl)*Vvec(il + 1);
            end 
        end

        function kgrid = getKgrid(nk, kl, kh)
            %make logpace kgrid
        
            kgrid = logspace(log(kl - kl +1)/log(10.0), log(kh - kl + 1)/log(10.0), nk)';
            kgrid = kgrid + ones(size(kgrid))*(kl - 1); 
            kgrid = kgrid';
        
        end

        function grid = logspace(l, h, nx)
            grid = logspace(log(l + -1.0*l + 1.0)/log(10.0), ...
                log(h + -1.0*l + 1.0)/log(10.0), nx)';
            grid = grid + ones(size(grid))*(l - 1.0);
        end

        function [P_mat, z_grid] = getTauchen(Nz, mu, sigma, rho, s)
            
            sigma_x = ((sigma^2)/(1-rho^2))^(0.5); 

            x_1 = mu - s*sigma_x;
            x_Nz = mu + s*sigma_x;
            
            x_grid = linspace(x_1, x_Nz,Nz);
            z_grid = exp(x_grid);
            
            %define interval width internally
            w = x_grid(1,2) - x_grid(1,1);
            
            %our transition matrix is Nz by Nz, with an algorithmic way of constructing
            %probabilities
            
            P_mat = zeros(Nz,Nz);
            
            for r = 1:Nz
                x_curr = x_grid(r)*rho;

                %fill in the first and last entries
                P_mat(r,1) = normcdf(x_grid(1) - x_curr + w/2, mu, sigma);
                P_mat(r,Nz) = 1 - normcdf(x_grid(Nz) - x_curr - w/2, mu, sigma);
%                 fprintf("Row %i: F1 = %1.4f; mean = %2.4f; std = %2.4f\n", ...
%                     r, F1, mu, sigma)

                for c = 2:Nz-1
                    upper = normcdf(x_grid(1,c) - x_curr + w/2, mu, sigma);
                    lower = normcdf(x_grid(1,c) - x_curr - w/2, mu, sigma);
                    P_mat(r,c) = upper - lower;
                end
            end
        end

        function [ix,we] = weight(pim, f)
            if f >= max(pim)
                ix = max(size(pim)) - 1;
                we = 0;
            elseif f < min(pim)
                ix = 1;
                we = 1;
            else
                ix = sum(pim <= f);
                we = (pim(ix+1) - f)/(pim(ix+1) - pim(ix));
            end
        end

        function [C] = cubicSpline(l, h, r, VK)
            lb = log(l)/log(10);
            ub = log(h)/log(10);
            
            val_grid = logspace(lb,ub,r+2);
            
            ffunc = @(ki) VK(1,ki);

            %make coeff matrix:
            
            ti = val_grid;
            fti = zeros(r+2,1);
            dti = zeros(r+2,1);
            ttf = zeros(r+2,1);
            
            fti(1) = ffunc(1);
            
            for i = 2:r+2
                t = val_grid(i);
                t_1 = val_grid(i-1);
            
                dti(i,1) = (t - t_1);
                fti(i,1) = ffunc(i);
                ttf(i,1) = (fti(i,1) - fti(i-1,1))/dti(i,1);
            end
            
            upper_diag = zeros(r+1,1);
            lower_diag = zeros(r+1,1);
            principal = zeros(r+1,1);
            
            for i = 2:r
                upper_diag(i,1) = dti(i,1);    
                lower_diag(i,1) = dti(i+2,1);
            end
            
            upper_diag(2,1) = upper_diag(2,1) + (dti(2,1)^2)/dti(3,1);
            lower_diag(r,1) = lower_diag(r,1) + (dti(r,1)^2)/dti(r-1,1);
            
            omega_1 = dti(3,1) - (dti(2,1)^2)/dti(3,1);
            omega_r = dti(r+1,1) - (dti(r+2,1)^2)/dti(r+1,1);
            
            principal(1,1) = 2*(dti(2,1) + dti(3,1)) - omega_1;
            principal(r+1,1) = 2*(dti(r+1,1) + dti(r+2,1)) - omega_r;
            
            for i = 2:r
                principal(i,1) = 2*(dti(i+1,1) + dti(i+2,1));
            end
            
            f = zeros(r,1);
            f(1,1) = 3*(dti(3,1)*(ttf(2,1)) + dti(2,1)*ttf(3,1)) - ...
                2*(dti(3)*ttf(2) - (dti(2)^2)/dti(3)*ttf(3,1));
            
            f(r,1) = 3*(dti(r+2)*(ttf(r+1)) + dti(r+1)*ttf(r+2)) - ...
                2*(dti(r+1)*ttf(r+2) - (dti(r+2)^2)/dti(r+1)*ttf(r+1));
            
            for i = 2:r
                f(i,1) = 3*(dti(i+2,1)*ttf(i+1,1) + dti(i+1)*ttf(i+2));
            end
            
            T = full(spdiags([circshift(lower_diag,-1) principal upper_diag], -1:1, r,r));
            T(r,r) = 2*(dti(r+1,1) + dti(r+2,1)) - omega_r;
            
            s = T\f;
            
            s_0 = 2*ttf(2,1)-((dti(2)/dti(3))^2)*ttf(3,1) - ...
                (1-(dti(2,1)/dti(3,1)^2))*s(1) + ...
                (dti(2,1)/dti(3,1))^2*s(2);
            
            s_r1 = 2*(ttf(r+2,1)-((dti(r+2,1)/dti(r+1,1)))^2*ttf(r+1)) - ...
                (1-(dti(r+2,1)/dti(r+1,1))^2)*s(r) + ...
                (dti(r+2,1)/dti(r+1,1))^2*s(r-1);
            
            s_fin = [s_0; s; s_r1];
            
            %finally, getting our coef matrix
            
            C = zeros(r+1, 4);
            
            for i = 1:r+1
                C(i,1) = fti(i);
                C(i,2) = s_fin(i,1);
                C(i,3) = 3*ttf(i+1)/dti(i+1) - ...
                    2*s_fin(i,1)/dti(i+1) - s_fin(i+1,1)/dti(i+1);
                C(i,4) = (-2*ttf(i+1) + ...
                    s_fin(i) + s_fin(i+1))/(dti(i+1)^2);
            end

        end

        function val = getSplineVal(coeffs, lval, lgrid)
            il = sum(lgrid <= lval);
            x = lgrid(il);
            if il == length(lgrid)
                il = il - 1;
            end

            d = lval - x;

            coeffs = squeeze(coeffs);

            val = coeffs(1,il) + ...
                coeffs(2, il)*(d) + ...
                coeffs(3,il)*(d)^2 + ...
                coeffs(4,il)*(d)^3;
        end
         
        function x = dist(M, N, nd)
            x = abs(M-N);
            for i = 1:nd
                x = max(x);
            end
        end

        function distr = condense(adistr, amu, agrid)
            [nl, nmu] = size(adistr);
            [na] = length(agrid);
            distr = zeros(nl, na);
            for im = 1:nmu
                [ix, we] = compute.weight(agrid, amu(im));

                distr(:, ix) = distr(:, ix) + we*adistr(:, im);
                distr(:, ix+1) = distr(:, ix+1) + (1-we)*adistr(:, im);
            end
        end
   
        function [TV, G, EV] = interpV(terms, V, EV, vTol)
            beta = terms.beta;
            sigma = terms.sigma;
            pil = terms.pil;
            agrid = terms.agrid;
            lgrid = terms.lgrid;
            na = length(agrid); nl = length(lgrid);
            r = terms.r;
            w = terms.w;
            lambda = terms.lamval;
            tau = terms.tau;
            captax = terms.captax;
            G = zeros(size(V)); TV = G;

            % now gss over spline
            for il = 1:nl
                l = lgrid(il);
                tau_r = captax(il);
                for ia = 1:na
                    a = agrid(ia);
                    y = gov.tax(w*l,lambda, tau) + (1+r*(1-tau_r))*a;
                    ix = sum(agrid <= y);
                    ix = max(ix, 1);  % ensure at least one index
                    searchgrid = agrid(1:ix);
                    [vval, aval] = compute.gss(EV(il,:),y, beta, ...
                        sigma, searchgrid, vTol*1e-2);
                    TV(il, ia) = vval;
                    G(il, ia) = aval;
                end
            end

            for ia = 1:na
                for il = 1:nl
                    EV(il, ia) = pil(il,:)*TV(:,ia);
                end
            end

        end

        function [V, G] = backsolve(terms,EV,vTol)

            beta = terms.beta;
            sigma = terms.sigma;
            pil = terms.pil;
            agrid = terms.agrid;
            lgrid = terms.lgrid;
            na = length(agrid); nl = length(lgrid);
            r = terms.r;
            w = terms.w;
            lambda = terms.lamval;
            tau = terms.tau;
            captax = terms.captax;
            G = zeros(size(EV)); V = G;

            % now gss over spline to get new V, G.
            for il = 1:nl
                l = lgrid(il);
                tau_r = captax(il);
                for ia = 1:na
                    a = agrid(ia);
                    y = gov.tax(w*l,lambda, tau) + (1+r*(1-tau_r))*a;
                    ix = sum(agrid <= y);
                    ix = max(ix, 1);  % ensure at least one index
                    searchgrid = agrid(1:ix);
                    [vval, aval] = compute.gss(EV(il,:),y, beta, ...
                        sigma, searchgrid, vTol*1e-2);
                    V(il, ia) = vval;
                    G(il, ia) = aval;
                end
            end

        end

    end

end