classdef aubhik
    methods(Static)

        function [c] = spline(r, dim, knots, fvals)
            
            % Set the dimension of the problem, the current example is univariate
            m = dim;
                                    
            % endpoint condition choices are either 'not-a-knot' or 'complete'
            % the default is 'not-a-knot' while the complete choice requires 
            % slope values for s(0) and s(r+1). SVEC contains these values 
            % which are irrelevant when the not-a-knot condition is chosen.
            SVEC = [0 0]';
            % indicator = 'complete';
            indicator = 'not-a-knot';
            
            % Compute spline arrays which are determined by knots
            % Note these matrices are independent of functions values
            [L, U, dtau] = aubhik.spa(knots, r, indicator);
            
            % Compute piecewise polynomial cubic spline interpolant
            [c] = aubhik.spb(knots, fvals, r, m, L, U, dtau, indicator, SVEC);

        end

        function [L, U, dtau] = spa(knots, r, indicator)

            % Aubhik Khan and Julia Thomas
            %
            % knots is a vector of length r+2 containing the knots
            % r is the number of interior knots points
            % indicator is a character string which specifies the endpoint condition
            %
            % This program sets up the LHS for cubic spline interpolation 
            % based on the not-a-knot condition.  This LHS is invariant to
            % changes in function values that may arise, for example, during
            % repeated applications of a contraction map that updates a 
            % value function. 
            %
            % program family for univariate spline calibration (creation)
            % spa (this program) (external input for L, U and dtau)
            % spb containing SPRHS, SPLSolve, SPUSolve, SPppcoeff (internal subroutines)
            % associated program for spline evaluation speval
            
            if (nargin < 3)
               indicator = 'not-a-knot';
            end
            
            % starting index (0 in notes) should be read as zero.
            z = 1;
            dtau = knots(z+1:r+z+1) - knots(z:r+z);
            
            % T is sparse matrix storage of the tridiagonal matrix T
            % T(1,:) is the upper diagonal, T(2,:) is the principal diagonal
            % and T(3,:) is the lower diagonal.  In other words, we store a 
            % r x r matrix as a 3 x r.  In this transformation the stored T(l,j)
            % corresponds to the true T(j-1,j), T(u,j) to T(j+1,j) and T(p,j) to
            % T(j,j).
            
            T = zeros(3,r);
            
            % principal, lower and upper diagonals
            p = 2; 
            l = 3;
            u = 1; 
            
            if (indicator(1:5) == 'not-a' )
               ratio0 = dtau(z)/dtau(z+1); 
               ratio0 = ratio0*dtau(z);
               nak11 = -dtau(z+1) + ratio0;
               nak12 = ratio0;
               
               ratior = dtau(z+r)/dtau(z+r-1); 
               ratior = ratior*dtau(z+r);
               nakrr = -dtau(z+r-1) + ratior;
               nakrrm1 = ratior;
            else
               nak11 = 	0;
               nak12 = 0;
               nakrr = 0;
               nakrrm1 = 0;
            end
               
            % initialize first row of T
            
            T(p,1) = nak11 + 2*(dtau(z) + dtau(z+1));
            T(u,1) = nak12 + dtau(z);
            
            % initialize first regular row
            
            T(l,2) = dtau(z+2);
            T(u,2) = dtau(z+1);
            T(p,2) = 2*(T(l,2) + T(u,2));
            
            % initialize remaining rows until final row
            
            for j = 3:1:r-1
               T(l,j) = dtau(z+j);
               T(u,j) = T(l,j-1);
               T(p,j) = 2*(T(l,j) + T(u,j));
            end
            
            % initialize final row
            
            T(p,r) = nakrr + 2*(dtau(z+r) + dtau(z+r-1));
            T(l,r) = nakrrm1 + dtau(z+r);
            
            % This is the actual matrix of coefficients that is being stored in sparse form as T.
            % TX = diag(T(p,:),0) + diag(T(u,1:r-1),1) + diag(T(l,2:r),-1);
            
            [L, U] = aubhik.TridiagLUfactor(T);
            
            % These are the actual matrices that represent the diagonalization of the tridiagonal matrix T.
            % LX = diag(L,-1) + diag(ones(1,r));
            % UX = diag(U(p,:),0) + diag(U(u,1:r-1),1);
            
            % sum(sum(abs(TX - LX*UX)))
        end 

        function [Lm,Um] = TridiagLUfactor(Tm)
        
                n = length(Tm);
                
                % The LU factorization is stored as a sparse matrix.  Only the 
                % lower diagonal of L is computed, the principal diagonal being 
                % ones. Both the principal and upper diagonal of Um are computed.
                
                % The stored Lm(j-1) is the true Lm(j,j-1).
                
                Lm = zeros(1,n-1); 
                Um = zeros(2,n);
                u = 1;
                p = 2;
                l = 3;
                Um(u,1:n-1) = Tm(u,1:n-1);
                
                Um(p,1) = Tm(p,1);
                
                for i = 2:1:n-1
                   Lm(i-1) = Tm(l,i)/Um(p,i-1);
                   Um(p,i) = Tm(p,i) - Lm(i-1)*Um(u,i-1);
                end
                
                Lm(n-1) = Tm(l,n)/Um(p,n-1);
                Um(p,n) = Tm(p,n) - Lm(n-1)*Um(u,n-1);
                
        end

        function [c] = spb(knots, fvals, r, m, L, U, dtau, indicator, SVEC)
        
            % Aubhik Khan and Julia Thomas
            %
            % knots is a vector of length r + 2 containing the knots
            % fvals is an array of dimension r x m, containing m function values at the knots
            % r is the number of interior knot points
            % m is the number of functions being interpolated over the same knots
            % L, U and dtau are function independent and given by a call to spa
            % indicator and SVEC specify the endpoint condition
            %
            % This program requires the LU factorization of the LHS to have been 
            % completed using spa.m and thus to have L, U and dtau as inputs.  
            % It also also requires knots, function values, the number of interior 
            % knots and the number of functions subject to spline approximation.
            %
            % program family for univariate spline calibration (creation)
            % spa (external input for L, U and dtau)
            % spb (this program) containing SPRHS, SPLSolve, SPUSolve, SPppcoeff (internal subroutines)
            % associated program for spline evaluation speval
            
            % Set up RHS in search for interior slopes.
            if (indicator(1:5) == 'not-a')
               [B, df] = aubhik.SPRHS(dtau', r, m, fvals');
            else
               [B, df] = aubhik.SPRHS(dtau', r, m, fvals', indicator, SVEC);
            end
            
            %An alternative to the factorization is simply s = TX\B;
            
            % Perform forward pass using L matrix to determine Ly = B
            [y] = aubhik.SPLsolve(L,B);
            
            % Perform backwards pass using U matrix to determine Ux = y
            [s] = aubhik.SPUsolve(U,y);
            
            % save Sppout fvals B df s y L U dtau 
            
            % s is 1 x r with the interior slopes.
            if (indicator(1:5) == 'not-a')
               [c] = aubhik.SPppcoeff(fvals, dtau, r, m, s, df);
            else
               [c] = aubhik.SPppcoeff(fvals, dtau, r, m, s, df, indicator, SVEC);
            end
        
        end

        function [f,df] = SPRHS(dtau, r, m, fvals, indicator, SVEC)
        
            % part of the SPpp and SPppmaster family of programs
            
            if (nargin < 5)
               indicator = 'not-a-knot';
            end
            
            % translate the zero index so that s(0) becbomes s(z).
            z = 1;
            
            %  compute divided difference of f: (f[tau(i+1)] - f[tau(i)])/dtau(i);
            df = fvals(z+1:z + r + 1,:) - fvals(z:z + r,:);
            df = df./dtau;
            
            f = zeros(r,m);
            
            if (indicator(1:5) == 'not-a')
               
               % These endpoint conditions are the implication of the not-a-knot condition
               ratio0 = dtau(z,:)./dtau(z+1,:);
               ratio0 = ratio0.*dtau(z,:);
               nak1 = -2*dtau(z+1,:).*df(z,:) + 2*ratio0.*df(z+1,:);
               
               ratior = dtau(r+z,:)./dtau(r+z-1,:);
               ratior = ratior.*dtau(r+z,:);
               nakr = 2*ratior.*df(r+z-1,:) - 2*dtau(r+z-1,:).*df(r+z,:);
               
            else
                  
               nak1 = -dtau(z+1,:).*SVEC(1,:);
               nakr = -dtau(z+r-1,:).*SVEC(2,:);
               
            end
            
            % Set first and last RHS terms given the endpoint conditions.
            f(1,:) = dtau(z+1,:).*df(z,:) + dtau(z,:).*df(z+1,:);
            f(1,:) = 3*f(1,:) + nak1;
            f(r,:) = dtau(r+z,:).*df(r+z-1,:) + dtau(r+z-1,:).*df(r+z,:);
            f(r,:) = 3*f(r,:) + nakr;
            
            for j = 2:1:r-1
               f(j,:) = dtau(z+j,:).*df(z + j - 1,:) + dtau(z + j - 1,:).*df(z + j,:);
               f(j,:) = 3*f(j,:);
            end
        
        end

        function [y] = SPLsolve(L,B)
        
            % forward pass for diagonalized tridiagonal spline matrix
            
            [n,r] = size(B);
            
            % L is of length(n-1) and is the lower diagonal of a bi-diagonal 
            % matrix described by trueL(j,j) = 1, trueL(j-1,j) = L(j-1) for j>1.
            y = zeros(n,r);
            
            y(1,:) = B(1,:);
            
            for i = 2:1:n
               y(i,:) = B(i,:) - L(i-1)*y(i-1,:);
            end
        
        end

        function [x] = SPUsolve(U,y)
            
            % backward pass for diagnolized tridiagonal matrix
            
            [n,r] = size(y);
            p = 2;
            u = 1;
            
            % U is a 2 x n vector that describes a true matrix U where trueU(j,j) = U(p,j) and 
            % trueU(j,j+1) = U(u,j) for j < n.
            
            x = zeros(n,r);
            
            x(n,:) = y(n,:)/U(p,n);
            
            for i = n-1:-1:1
               x(i,:) = (y(i,:) - U(u,i)*x(i+1,:))/U(p,i);
            end
        
        end

        function [c] = SPppcoeff(fvals, dtau, r, m, s, df, indicator, SVEC)
        
            % This retrieves the c(1) - c(4) coefficients of the interpolating cubic spline given 
            % that the slopes s(1) - s(r) that have already been computed elsewhere.
            
            if (nargin < 7) 
               indicator = 'not-a-knot';
            end
            
            dtau2 = dtau.*dtau;
            
            % translate the zero starting index value
            z = 1;
            nc = (r+z)*m;
            c = zeros(4,nc);
            
            if (indicator(1:5) == 'not-a' )
               
               % recover s(0) and s(r+1) given s(1),...,s(r) 
               % for the not-a-knot endpoint conditions that 
               % imply third derivative continuity at the 
               % first and last interior knots.
               
               dratios0 = dtau(z,:)./dtau(z+1,:);
               dratios0 = dratios0.^2;
               dratiosrp1 = dtau(z+r,:)./dtau(z+r-1,:);
               dratiosrp1 = dratiosrp1.^2;
               
               s0 = (-1 + dratios0).*s(1,:);
               s0 = s0 + dratios0.*s(2,:) + 2*df(z,:) - 2*dratios0.*df(z+1,:);
               
               srp1 = (-1 + dratiosrp1).*s(r,:);
               srp1 = srp1 + dratiosrp1.*s(r-1,:) + 2*df(z+r,:) - 2*dratiosrp1.*df(z+r-1,:);
               
            else
               
               % assume the complete spline
               s0 = SVEC(1,:);
               srp1 = SVEC(2,:);
               
            end
            
            % augmented slope vector
            sadj = [s0; s; srp1];
            
            c(1,:) = reshape(fvals(z:r+z,:),1,nc);
            c(2,:) = reshape(sadj(z:r+z,:),1,nc);
            
            c(3,:) = reshape((df(z:r+z,:)./dtau(z:r+z,:)),1,nc);
            c(3,:) = 3.*c(3,:) - 2.*reshape(sadj(z:r+z,:)./dtau(z:r+z,:),1,nc);
            c(3,:) = c(3,:) - reshape(sadj(z+1:r+z+1,:)./dtau(z:r+z,:),1,nc);
            
            c(4,:) = -reshape(df(z:r+z,:)./dtau2(z:r+z,:),1,nc);
            c(4,:) = 2*c(4,:) + reshape(sadj(z:r+z,:)./dtau2(z:r+z,:),1,nc);
            c(4,:) = c(4,:) + reshape(sadj(z+1:r+z+1,:)./dtau2(z:r+z,:),1,nc);
        
        end

    end
end