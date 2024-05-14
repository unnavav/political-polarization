classdef cubic
    methods(Static)

        function [c] = spline(r, dim, knots, fvals)
            % Set the dimension of the problem, the current example is univariate
            m = dim;
                                    
            % Endpoint condition choices are either 'not-a-knot' or 'complete'
            indicator = 'not-a-knot';
            SVEC = [0 0]';
            
            % Compute spline arrays which are determined by knots
            [L, U, dtau] = cubic.spdiag(knots, indicator);
            
            % Compute piecewise polynomial cubic spline interpolant
            c = cubic.spcoeff(fvals, r, L, dtau', indicator, SVEC);
        end

        function [L, U, dtau] = spdiag(knots, indicator)
            % Function to set up the LHS for cubic spline interpolation based on endpoint condition
            dtau = diff(knots);

            % Create tri-diagonal matrix for cubic spline system
            L = diag(2 * (dtau(1:end-1) + dtau(2:end))) + diag(dtau(2:end-1), 1) + diag(dtau(2:end-1), -1);
            U = diag(dtau(2:end), 1);

            % Include endpoint conditions (assuming 'not-a-knot')
            if strcmp(indicator, 'not-a-knot')
                L(1, 1:2) = [1, -1];
                L(end, end-1:end) = [-1, 1];
            end
        end

        function [c] = spcoeff(fvals, r, L, dtau, indicator, SVEC)
            % Compute coefficients for piecewise polynomial spline
            df = diff(fvals) ./ dtau;
            s = L \ (6 * (df(2:end) - df(1:end-1)));
            
            % Include endpoint conditions (assuming 'not-a-knot')
            if strcmp(indicator, 'not-a-knot')
                s = [0; s; 0];
            else
                s = [SVEC(1); s; SVEC(2)];
            end
            
            % Compute coefficients for piecewise polynomials
            c = zeros(4, r + 1);
            c(1, :) = fvals(1:end-1);
            c(2, :) = s(1:end-1) .* dtau;
            c(3, :) = 3 * (fvals(2:end) - fvals(1:end-1)) ./ dtau - 2 * s(1:end-1) - s(2:end);
            c(4, :) = 2 * (fvals(1:end-1) - fvals(2:end)) ./ dtau + s(1:end-1) + s(2:end);
        end

    end
end
