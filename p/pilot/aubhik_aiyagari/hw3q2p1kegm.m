function [termsout] = hw3q2p1kegm(termsarray, r, wage, mu0, v0)

% mu, kagg, v, g, muGI
% termsout.mu = mu;
% termsout.kagg = kagg;
% termsout.v = v;
% termsout.g = g;
% termsout.muGI = muGI;


% 20 November 2020 I followed In Hwan and added distance2 = abs(tkagg -
% kagg)) in line 240.

% lines 59 - 66 are important to understand how egm handles borrowing
% limits.  26 October 2015.

%flnum = 'lnum'; flgrid = 'lgrid'; fpil = 'pil'; fbeta = 'beta'; fsigma = 'sigma'; fanum = 'anum'; fagrd = 'agrid';
%fphi = 'phi'; fmuagnum = 'muagnum'; fmuagrid = 'muagrid'; fprecision = 'precision'; fprecision1 = 'precision1';

lnum = termsarray.lnum; lgrid = termsarray.lgrid; pil =  termsarray.pil; beta = termsarray.beta; sigma = termsarray.sigma; 
anum = termsarray.anum; agrid = termsarray.agrid; phi = termsarray.phi; muagnum = termsarray.muagnum; muagrid = termsarray.muagrid;
precision = termsarray.precision; precision1 = termsarray.precision1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           %
%       Step 1: Solve the value function    %
%                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method of successive approximations
v = zeros(size(v0)); g = v; tv = v; tg = g;

ahategm = v;

iteration = 0; distance1 = 2.0*precision;

ev = zeros(lnum, anum);

v = v0;

disp ( ' ' )
fprintf( ' solving (r,w) = (%8.4f, %8.4f) with endogenous grid \n ', r, wage );

while (distance1 > precision )
    
    % conditional expectation v(e,b') = E{v(e',b')|e}
    for m = 1:1:lnum
        ev(m,1:anum) = pil(m,1:lnum)*v(1:lnum,1:anum);
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
    
    %disp ( ' ahategm ' )
    %disp(ahategm)
    %pause
    
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
            g(m,ia) = afval;
            
        end
        
    end
    
    % disp ( ' g ' )
    % disp(g)
    % pause
    
    % Update the value function using the decision rule
    for ia = 1:1:anum
        ahat = agrid(ia);
        for m = 1:1:lnum
            lval = lgrid(m);
            
            % The choice set for qb' must be bounded between blow and min(bhigh, yval).
            zval = wage*lval + (1.0 + r)*ahat - r*phi;
            afval = g(m,ia); cval = zval - afval;
            [objective] = operator(sigma, beta, anum, agrid, afval, ev(m, 1:anum), cval);
            tv(m,ia) = objective;
        end
    end   
    
  
    distance = max(max((abs(tv - v))));
    v = tv;
    distanceg = max(max(abs(tg - g)));
    tg = g;
    distance1 = min(distance, distanceg);
    
    amin = min(min(g)); amax = max(max(g));
    iteration = iteration + 1;
    
    if (mod(iteration,250) == 0) 
        fprintf ( ' iter %4d    ||Tv-v|| = %8.4f      ||tg - g|| = %8.4f     ahat (min, max) = (%8.4f, %8.4f) \n ', ...
            iteration, distance, distanceg, amin, amax);
    end
end

fprintf ( ' iter %4d    ||Tv-v|| = %8.4f      ||tg - g|| = %8.4f     ahat (min, max) = (%8.4f, %8.4f) \n ', ...
            iteration, distance, distanceg, amin, amax);
        
save valueaegm v g


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                           %
%       Step 3: Compute stationary distribution             %
%                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = zeros(lnum, muagnum); muGI = mu; muGIweight = mu;
solve = mu;

distance1 = 2.0*precision1; 
distance2 = distance1;
iteration = 0;

% 13 November 20, check covergence in kagg and mu after In Hwan's findings this summer.
mu = mu0;
kagg = 0;

%disp ( ' ' )
%fprintf ( ' solving for a stationary distribution \n ' )

while (max(distance1, distance2)  > precision1)

    muf = zeros(lnum, muagnum);

    for i = 1:1:muagnum
        ahat = muagrid(i);
        for m = 1:1:lnum
         
            muval = mu(m,i);
            if (muval > 0)

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                                                           %
                %     STEP 3.1: DECISION RULES ON GRID FOR DISTRIBUTION     %
                %                                                           %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if (solve(m,i) == 0)
                    
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
                    afhat = g(m,aindex)*weight + g(m,aindex+1)*(1.0 - weight);
                    
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
                    
                    %s = sprintf( ' ahat = %8.4f    afhat = %8.4f   aindex = %4d    muagrid(aindex) = %8.4f    muagrid(aindex+1) = %8.4f   weight = %8.4f   weighted value = %8.4f ', ...
                    %    ahat, afhat, aindex, muagrid(aindex), muagrid(aindex+1), weight, weight*muagrid(aindex) + (1.0 - weight)*muagrid(aindex+1));
                    %disp (s)
                    muGI(m,i) = aindex;
                    muGIweight(m,i) = weight;
                    solve(m,i) = 1;
                    
                else
                    aindex = muGI(m,i); weight = muGIweight(m,i);                                    
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

    tmu = muf;

    tkagg = 0.0;
    for m = 1:1:lnum
        tkagg = tkagg + mu(m,1:muagnum)*muagrid';
    end
    tkagg = tkagg - phi;

    % 13 November 20, make sure kagg has converged.
    distance1 = max(max(abs(tmu - mu)));
    distance2 = abs(tkagg - kagg);

    iteration = iteration + 1;
    
    if (mod(iteration,750) == 0)
        fprintf( ' iter %4d    (distmu, distk) = (%8.6f, %8.6f)   sum = %6.4f   tkagg = %8.4f \n ', iteration, distance1, distance2, sum(sum(mu)), tkagg);
     end

    mu = tmu;
    kagg = tkagg;
 

end

 fprintf( ' iter %4d    (distmu, distk) = (%8.6f, %8.6f)   sum = %6.4f   tkagg = %8.4f \n ', iteration, distance1, distance2, sum(sum(mu)), tkagg);

% structure to return results
% mu, kagg, v, g, muGI
termsout.mu = mu;
termsout.kagg = kagg;
termsout.v = v;
termsout.g = g;
termsout.muGI = muGI;


%%%%%%%%%%%%%%%%%%%%%%%%%
%                       %
%       Operator        %
%                       %
%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply the operator having found the decision rule using egm.
function [objective] = operator(sigma, beta, anum, agrid, afval, ev, cval)

aindex = gridlookup(anum, agrid, afval);

weight = agrid(aindex+1) - afval;
weight = weight/(agrid(aindex+1) - agrid(aindex));
weight = max(weight, 0);
weight = min(weight, 1);

evalf = ev(aindex)*weight + ev(aindex+1)*(1.0 - weight);

if (sigma == 1.0)
    fval = log(cval) + beta*evalf;
else
    fval = (cval^(1.0-sigma))/(1.0 - sigma) + beta*evalf;
end

objective = fval;
