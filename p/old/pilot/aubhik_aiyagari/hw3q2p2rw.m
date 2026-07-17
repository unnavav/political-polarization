function [kaggvec] = hw3q2p2rw(termsarray, tmax, rvec, wagevec, mu0, vss)

% adapted from aiyagkegm Aubhik 5 March 2023
% this solves for a perfect foresight path given rvec(t) and wagevec(t)
%
% it's very fast to write from aiygagkegm which solves for v and mu in the
% steady state.  This program adapts it to just iterate over time, and for
% r and wage that vary with t, and not until convergence.  

% lnum = termsarray.lnum; lgrid = termsarray.lgrid; pil =  termsarray.pil; beta = termsarray.beta; sigma = termsarray.sigma; 
% anum = termsarray.anum; agrid = termsarray.agrid; phi = termsarray.phi; muagnum = termsarray.muagnum; muagrid = termsarray.muagrid;
% precision = termsarray.precision; precision1 = termsarray.precision1; 
lnum = termsarray.lnum;
lgrid = termsarray.lgrid;
pil = termsarray.pil;
beta = termsarray.beta;
sigma = termsarray.sigma;
anum = termsarray.anum;
agrid = termsarray.agrid;
phi = termsarray.phi;
muagnum = termsarray.muagnum;
muagrid = termsarray.muagrid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           %
%       Step 1: Solve the value function    %
%                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% value function for each period, tmax is steady state.
v = zeros(lnum, anum, tmax);
g = v;

v(:,:,tmax) = vss;

ahategm = zeros(lnum, anum);
ev = zeros(lnum, anum);

fprintf ( ' solving value functions with endogenous grid \n ');
for t = tmax-1:-1:1
    
    r = rvec(t); wage = wagevec(t);
    tv = v(:,:,t+1);

    % conditional expectation v(e,b') = E{v(e',b')|e}
    for m = 1:1:lnum
        ev(m,1:anum) = pil(m,1:lnum)*tv(1:lnum,1:anum);
    end

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
    
    for iaf = 1:1:anum
        afhat = agrid(iaf);
        
        for m = 1:1:lnum
            lval = lgrid(m);
            Devaf = Dev(m,iaf);

            zval = wage*lval - r*phi;
       
            [ahat] = decisionegm(sigma, beta, afhat, Devaf, r, zval); 
            ahategm(m,iaf) = ahat;
           
        end
    end
    
    % endogenous grid method piecewise linear interapolation of decision
    % rules forward, afhat = g(lval, ahat).
    for m = 1:1:lnum
        aegmgrid(1:anum) = ahategm(m,1:anum);
        for ia = 1:1:anum
            ahat = agrid(ia);
            
            % We know that afhat grid began with agrid(1).  If the first 
            % endogenous grid point aegmgrid(1) exceeds ahat, as
            % agrid(1) is the decision for aegmgrid(1), monotonicity of the 
            % decision rule given the borrowing limit says that it is also  
            % the decision for ahat.  This is how the borrowing limit 
            % is implemented.
            if (ahat < aegmgrid(1)) 
                afval = agrid(1);
            else
            
                aindex = gridlookup(anum, aegmgrid, ahat);

                if(aindex < anum)
                    weight = aegmgrid(aindex+1) - ahat;
                    weight = weight/(aegmgrid(aindex+1) - aegmgrid(aindex));
                    weight = max(weight, 0);
                    weight = min(weight, 1.0);
                else
                    aindex = aindex - 1;
                    weight = 0.0;
                end

                % monotonicity of decision rule
                afval = weight*agrid(aindex) + (1.0 - weight)*agrid(aindex + 1);

                % borrowing limit is simply imposed in endogenous grid method
                % approach using monotonicity of decision rule.
                if (afval < agrid(1))
                    afval = agrid(1);
                end
            
            end
            g(m,ia,t) = afval;
            
        end
        
    end
    
    
    % Update the value function using the decision rule
    for ia = 1:1:anum
        ahat = agrid(ia);
        for m = 1:1:lnum
            lval = lgrid(m);
            
            % The choice set for qb' must be bounded between blow and min(bhigh, yval).
            zval = wage*lval + (1.0 + r)*ahat - r*phi;
            afval = g(m,ia,t); cval = zval - afval;
            [objective] = operator(sigma, beta, anum, agrid, afval, ev(m, 1:anum), cval);
            v(m,ia,t) = objective;
        end
    end   
    
  
    distance = max(max((abs(tv - v(:,:,t)))));
    
    amin = min(min(g(:,:,t))); amax = max(max(g(:,:,t)));
    
    if (mod(t,50) == 0) 
        fprintf ( ' t = %4d    ||vf-v|| = %8.4f      ahat (min, max) = (%6.2f, %6.2f) \n ', ...
            t, distance, amin, amax);
    end
end
fprintf ( ' t = %4d    ||vf-v|| = %8.4f      ahat (min, max) = (%6.2f, %6.2f) \n ', ...
    t, distance, amin, amax);
fprintf ( ' \n ')

save valueimp v g


kaggvec = zeros(1,tmax-1);

muvec = zeros(lnum, muagnum, tmax);
muvec(:,:,1) = mu0;

mu = muvec(:,:,1);

for t = 1:1:tmax - 1

    for m = 1:1:lnum
        kaggvec(t) = kaggvec(t) + mu(m,1:muagnum)*muagrid';
    end
    kaggvec(t) = kaggvec(t) - phi;

    muf = zeros(lnum, muagnum);

    for ia = 1:1:muagnum
        ahat = muagrid(ia);

        for m = 1:1:lnum
         
            muval = mu(m,ia);
            if (muval > 0)

                    % We linearly interpolate the policy function on agrid to compute
                    % afhat on the distribution support.
                    
                    aindex = gridlookup(anum, agrid, ahat);
                    if(aindex < anum)
                        weight = agrid(aindex+1) - ahat;
                        weight = weight/(agrid(aindex+1) - agrid(aindex));
                        weight = max(weight, 0.0);
                        weight = min(weight, 1.0);
                    else
                        aindex = aindex - 1;
                        weight = 0.0;
                    end
                    afhat = g(m,aindex,t)*weight + g(m,aindex+1,t)*(1.0 - weight);
                    
                    aindex = gridlookup(muagnum, muagrid, afhat);
                    
                    if(aindex < muagnum)
                        weight = muagrid(aindex+1) - afhat;
                        weight = weight/(muagrid(aindex+1) - muagrid(aindex));
                        weight = max(weight, 0.0);
                        weight = min(weight, 1.0);
                    else
                        aindex = aindex - 1;
                        weight = 0.0;
                    end

                for mf = 1:1:lnum
                    if (aindex < muagnum)
                        muf(mf,aindex) = muf(mf,aindex) + pil(m,mf)*muval*weight;
                        muf(mf, aindex+1) = muf(mf, aindex+1) + pil(m,mf)*muval*(1.0 - weight);
                    else
                        muf(mf,aindex) = muf(mf,aindex) + pil(m,mf)*muval;
                    end
                end
               
            end

        end
    end

    muvec(:,:,t+1) = muf;

    mu = muf;
    
    if (mod(t,50) == 0 || t == 1)
        fprintf( ' t = %4d    kagg = %8.4f \n ', t, kaggvec(t));
    end
 
end

fprintf( ' t = %4d    kagg = %8.4f \n ', t, kaggvec(t));

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
%       subfunction for endogenous grid method      %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ahat] = decisionegm(sigma, beta, afval, Dev, r, zval)


% find consumption using inverse utility function and marginal expected
% value.
cval = (beta*Dev)^(-1.0/sigma);

ahat = (cval - zval + afval)/(1 + r);

%s = sprintf ( ' Dev = %8.4f  c = %8.4f  a = %8.4f ', Dev, cval, ahat);
%disp (s)

%%%%%%%%%%%%%%%%%%%%%%%%%
%                       %
%       Operator        %
%                       %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply the operator having found the decision rule using egm.
function [objective] = operator(sigma, beta, anum, agrid, afval, ev, cval)

aindex = gridlookup(anum, agrid, afval);

if(aindex < anum)
    weight = agrid(aindex+1) - afval;
    weight = weight/(agrid(aindex+1) - agrid(aindex));
    weight = max(weight, 0);
    weight = min(weight, 1);
else
    weight = 0.0;
    aindex = aindex - 1;
end

evalf = ev(aindex)*weight + ev(aindex+1)*(1.0 - weight);

if (sigma == 1.0)
    fval = log(cval) + beta*evalf;
else
    fval = (cval^(1.0-sigma))/(1.0 - sigma) + beta*evalf;
end

objective = fval;
