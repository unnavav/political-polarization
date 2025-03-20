classdef spline
    methods(Static)

        % Full spline solve--using interp1 for now, but will replace with
        % my own spline solution once I'm sure what I'm doing is actually
        % stable
        function [TV, G] = solve(terms, V, EV, vTol)
            beta = terms.beta;
            sigma = terms.sigma;
            pil = terms.pil;
            agrid = terms.agrid;
            lgrid = terms.lgrid;
            na = length(agrid); nl = length(lgrid);
            r = terms.r;
            w = terms.w;
            lambda = terms.lambda;
            tau = terms.tau;
            captax = terms.captax;
            aspline = compute.logspace(agrid(1), agrid(na), na*2);

            G = zeros(size(V)); TV = G;

            % now fit spline to expected value
            Vspline = zeros(nl, length(aspline));
            for il = 1:nl
                Vspline(il, :) = interp1(agrid,EV(il,:),aspline,'pchip')';
            end

            % now gss over spline
            for il = 1:nl
                l = lgrid(il);
                tau_r = captax(il);
                for ia = 1:na
                    a = agrid(ia);
                    y = w*l - gov.tax(w*l,lambda, tau) + (1+r*(1-tau_r))*a;
                    searchgrid = aspline(1:sum(aspline<=a));
                    [vval, aval] = compute.gss(Vspline(il,:),y, beta, ...
                        sigma, searchgrid, vTol*1e-2);
                    TV(il, ia) = vval;
                    G(il, ia) = aval;
                end
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
    end
end