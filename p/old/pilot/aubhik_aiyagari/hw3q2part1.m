% The Aiyagari (1994) model 
% S. Rao Aiyagari (1994) 'Uninisured idiosyncratic risk and aggregate
% savings' Quarterly Journal of Economics (August): 659 - 684.
%
% Aubhik Khan 2023
%
% Adapted Aiyagariphi1 for endogenous grid method.

clear
close all

% See the comments in Aiyagari (1994), pages 668- 669 and 671.  Equilibrium
% must have r* = 1/beta - 1 if there is no earnings risk.  Further, with
% idiosyncratic risk, kagg will go to infinity with infinitely lived
% households if r.ge.r*

% If ahigh is not sufficiently large, then the model will not solve for phi
% positive.  As rlow > rhigh, flow must be positive and fhigh negative for 
% bisection.  To solve phi = 2, ahigh must be 5.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%       Model parameters        %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A period is one year.  

% parameters
beta = 0.96173;         % subjective discount factor
sigma = 1.0;            % inverse of the elasticity of intertemporal substitution

alpha = 0.36;
delta = 0.08;
atfp = 1.0;

phi = 1.0;               % the borrowing limit is -phi.

% grid on labour productivity
lnum = 7;

% my calculations see quarter_annual
rhol = 0.9554;
stdl = 0.5327;

precision = 1.0e-6;     % numerical precision parameter for value function
factordist = 10;        % In Hwan remarks increased this value.  
                        % numerical precision parameter for distribution
precision1 = precision/factordist;

% discretised stochastic process for labour productivity
if (lnum > 1 )

    stdinnovl = stdl*sqrt((1.0 - rhol^2)); meanl= 0; multiple = 2.575;
    [lgrid, pil] = tauchen2(meanl, stdinnovl, rhol, multiple, lnum);
    lgrid = exp(lgrid);

    % aggregate efficiency units of labour
    pilr = ones(1, lnum)*(1/lnum);

    %vaasavi : iterating to convergence of the stationary distribution
    distance = 2.0*precision;
    while (distance > precision/1000.0)

        tpilr = (pilr*pil);
        distance = max(abs(tpilr - pilr));
        pilr = tpilr;
    end

else

    lgrid = 1;
    pil = 1;
    pilr = 1;
end

lagg = pilr*lgrid';

aphi = num2str(phi);


resultfile = 'resegmag';
resultfile = [resultfile 'phi' aphi];

% new way to time program
starttime = datetime('now');

% complete markets real interest rate
rstar = 1.0/beta - 1.0;

% initial bounds for real interest rate between 0 and 1/beta - 1
% implies a minimum capital stock that may exist in equilibrium.  
klow0 = (rstar + delta)/(alpha*atfp);
klow0 = klow0/(lagg^(1.0 - alpha));
klow0 = klow0^(1.0/(alpha - 1.0));

% The model is highly nonlinear in klow and khigh.  klow must be greater
% than the complete markets value of klow0 associated with rstar.  Next
% khigh must be very close.  At the same time there must be enormous
% precision in agrid and muagrid.
kl = 1.3;
kh = 1.5;

if (phi == 1)
    
    kl = 1.025;
    kh = 1.2;
    
elseif (phi >= 3)
    
    kl = 1.1;     % 13 November 20 must be low.
    kh = 1.4;
    
    if (rhol < 0.9)
        kl = 1.005;
        kh = 1.03;
    end
    
end

klow = kl*klow0;
khigh = klow0*kh;

% grid on modified assets, see Aiyagari (1994) equations 3a, 3b, 4a and 4a
% ahat = a + phi.  agrid is defined on ahat, hence its lower bound is 0.  
% This ensures the borrowing constaint is satisfied for any phi.  
% z = w*l + (1+r)*ahat - r*phi  and c = z - afhat.
alow = 0.0; ahigh = 100.0; anum = 250;
agrid = logspace(log(alow + -1.0*alow + 1.0)/log(10.0), log(ahigh + -1.0*alow + 1.0)/log(10.0), anum)';
agrid = agrid + ones(size(agrid))*(alow - 1.0);

%agrid = linspace(alow, ahigh, anum)';

% agrid = linspace(alow, ahigh, anum);
muagnum = 2500;
muagrid = linspace(alow, ahigh, muagnum);
%muagnum = anum; muagrid = agrid';

% phi = 3, ahigh = 40, anum = 100, muagnum = 2500, klow = 1.02 and khigh =
% 1.06

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%   Preliminary computation     %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% initial value function
v0 = zeros(lnum, anum); 

wage = 1.1; r = 0.04;

for i = 1:1:anum
    ahat = agrid(i);
    for m = 1:1:lnum
        lval = lgrid(m);
        % The choice set for qb' must be bounded between blow and min(bhigh, yval).
        zval = wage*lval + (1 + r)*ahat - r*phi;
        
        if (sigma == 1.0)
            v0(m,i) = log(zval)/(1.0 - beta);
        else

            v0(m,i) = (zval^(1 - sigma))/((1.0 - beta)*(1.0 - sigma));
        end

    end
end

disp ( ' % The model is discontinuous in klow and khigh.  klow must ' )
disp ( ' % be greater than the complete markets value of klow0 ' )
disp ( ' % associated with rstar. ' )

disp ( ' ' )
disp ( ' ** initial value function has to be increasing in l and a ** ' )
disp ( ' ' )

if (exist('valueaegm.mat', 'file')==2) 
    askme = input( ' use existing value function ? (y/n)', 's');
    if (askme == 'y') 
        load valueaegm
        if (size(v) == size(v0))
            v0 = v;
        end
    end
    
else
    disp ( ' press entry to continue  ' )
    pause
end


% simulate the Markov Chain
tlast = 10000;
cumpil = cumsum(pil, 2);
isimvec = zeros(1, tlast);
randomvariable = rand(1, tlast);

isimvec(1) = floor((lnum + 1)/2);
unit = ones(1, lnum);

for t = 2:1:tlast
    [no, isimvec(t)] = min( abs(cumpil(isimvec(t-1), 1:lnum) - randomvariable(t)));
end

lsimvec = lgrid(isimvec);

stdlsim = std(lsimvec);
rholsim = corrcoef(lsimvec(2:tlast), lsimvec(1:tlast-1)); rholsim = rholsim(1,2);

fprintf( ' sigma = %8.4f   stdl = %6.4f   rhol = %6.4f   simulated Markov process stdsiml = %6.3f   rhosiml = %6.3f  \n', sigma, stdl, rhol, stdlsim, rholsim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%       equilibrium real interest rate bisection        %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flnum = 'lnum'; flgrid = 'lgrid'; fpil = 'pil'; fbeta = 'beta'; fsigma = 'sigma'; fanum = 'anum'; fagrid = 'agrid';
fphi = 'phi'; fmuagnum = 'muagnum'; fmuagrid = 'muagrid'; fprecision = 'precision'; fprecision1 = 'precision1';

termsarray = struct(flnum, lnum, flgrid, lgrid, fpil, pil, fbeta, beta, fsigma, sigma, fanum, anum, fagrid, agrid, ...
    fphi, phi, fmuagnum, muagnum, fmuagrid, muagrid, fprecision, precision, fprecision1, precision1);


iteration = 0;
distance = 2*precision;
mu0 = zeros(lnum, muagnum);
mu0(1:lnum, floor(muagnum/3)) = pilr';

kval = klow;

r = alpha*atfp*(kval^(alpha - 1.0))*(lagg^(1.0 - alpha)) - delta;
wage = (1.0 - alpha)*atfp*(kval^alpha)*(lagg^(-alpha));
[termsout] = hw3q2p1kegm(termsarray, r, wage, mu0, v0);

% termsout structure with mu, kagg, v, g, and muGI
%mu = termsout.mu;
kagg = termsout.kagg;
v = termsout.v;
%g = termsout.g;
%muGI = termsout.muGI;

flow = kagg - kval;

fprintf( ' %4d     rstar = %8.4f  r = %8.4f   klow = %8.4f   kagg = %8.4f   flow = %8.4f \n', iteration, rstar, r, kval, kagg, flow);

if (kagg < 0) 
    fprintf( ' klow = %8.4f  r = %8.4f  kagg = %8.4f  \n', klow, r, kagg);
    pause 
end

v0 = v;
mu0(1:lnum, floor(muagnum/3)) = pilr';

kval = khigh;

r = alpha*atfp*(kval^(alpha - 1.0))*(lagg^(1.0 - alpha)) - delta;
wage = (1.0 - alpha)*atfp*(kval^alpha)*(lagg^(-alpha));
[termsout] = hw3q2p1kegm(termsarray, r, wage, mu0, v0);

% termsout structure with mu, kagg, v, g, and muGI
kagg = termsout.kagg;
v = termsout.v;

fhigh = kagg - kval;

fprintf( ' %4d     rstar = %8.4f  r = %8.4f  khigh = %8.4f   kagg = %8.4f  fhigh = %8.4f \n', iteration, rstar, r, kval, kagg, fhigh);

if (kagg < 0) 
    fprintf( ' khigh = %8.4f  r = %8.4f  kagg = %8.4f  \n', khigh, r, kagg);
    pause 
end

v0 = v;

if (flow*fhigh > 0.0)
    fprintf ( ' could not bisect for equilibrium real interest rate between (klow, khigh) = (%6.3f, %6.3f)   (kagg0, kagg1) = (%8.4f, %8.4f) \n', klow, khigh, flow + klow, fhigh + khigh);
else
    
    while (distance > precision/100.0)

        kval = (klow + khigh)/2.0;
        mu0 = zeros(lnum, muagnum);
        mu0(1:lnum, floor(muagnum/3)) = pilr';
        v0 = v;
        
        r = alpha*atfp*(kval^(alpha - 1.0))*(lagg^(1.0 - alpha)) - delta;
        wage = (1.0 - alpha)*atfp*(kval^alpha)*(lagg^(-alpha));
        [termsout] = hw3q2p1kegm(termsarray, r, wage, mu0, v0);

        % termsout structure with mu, kagg, v, g, and muGI
        kagg = termsout.kagg;
        v = termsout.v;

        fval = kagg - kval;
        
        if (kagg < 0)
            fprintf( ' kval = %8.4f  r = %8.4f  kagg = %8.4f  \n', kval, r, kagg);
            pause
        end

        if (fval*flow > 0.0)
            klow = kval; flow = fval;
        else
            khigh = kval; fhigh = fval;
        end

        distance  = khigh - klow;
        iteration = iteration + 1;

        fprintf( ' %4d     rstar = %8.4f   r = %8.4f  kval = %8.4f   kagg = %8.4f   fval = %8.4f \n', iteration, rstar, r, kval, kagg, fval);
                
    end
end

mu = termsout.mu;
g = termsout.g;
muGI = termsout.muGI;

% Marginal Propensity to Consume

% Per Krusell seminar 18 November said that the mpc for a representative
% agent model is dc/da when c = [r - delta]*a + w*n as a purely transitory
% shock is a rise in wealth.  

mpc = zeros(muagnum, lnum);
cons = zeros(muagnum, lnum);

for i = 1:1:muagnum
    
    % this is the wealth level from grid for the distribution.  
    ahat = muagrid(i);
    
    % compute grid lookup for ahat only if some labour productivity level
    % has positive number of people at ahat.  
    if (max(mu(1:lnum,i)) > 0)
        
        % look up ahat from muagrid on the agrid for the decision rules.
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
        
        % See the wealth level just above ahat, to eventually find the 
        % consmption change from a transitory shock which is captured as 
        % a rise in wealth.  
        if (i < muagnum)
            
            ahat1 = muagrid(i+1);  
            aindex1 = gridlookup(anum, agrid, ahat1);
            
            if(aindex1 < anum)
                weight1 = agrid(aindex1+1) - ahat1;
                weight1 = weight1/(agrid(aindex1+1) - agrid(aindex1));
                weight1 = max(weight1, 0.0);
                weight1 = min(weight1, 1.0);
            else
                aindex1 = aindex1 - 1;
                weight1 = 0.0;
            end            
        else            
            aindex1 = aindex;
            weight1 = weight;            
        end
        
    end
    
    % labour productivity levels
    for m = 1:1:lnum
        
        lval = lgrid(m);
        
        % number of people at (lval, ahat).
        muval = mu(m,i);
        
        if (muval > 0)
            % zval is cash on hand, remember that ahat is a - phi.  
            afhat = g(m,aindex)*weight + g(m,aindex+1)*(1.0 - weight);
            zval = wage*lval + (1.0 + r)*ahat - r*phi;
            cval = zval - afhat;
            
            afhat1 = g(m,aindex1)*weight1 + g(m,aindex1+1)*(1.0 - weight1);
            zval1 = wage*lval + (1.0 + r)*ahat1 - r*phi;
            cval1 = zval1 - afhat1;
            
            % marginal propensity to consume is the finite difference
            % approximation to the derivative of cval with respect to ahat.
            %  
            cons(i,m) = cval;
            
            if (ahat1 > ahat)
                mpc(i,m) = (cval1 - cval)/(ahat1 - ahat);
            else
                mpc(i,m) = mpc(i-1,m);
            end
        end
    end
    
end

meanmpc = sum(sum(mpc.*mu'));

% one dimensional mpc and mu over all lval.
mpcall = sum(mpc.*mu',2);
muall = sum(mu',2);
   
realrate = r;
realrate = realrate*100;

[Lmatrix, Amatrix] = meshgrid(lgrid, agrid - phi*ones(size(agrid)));

figure

subplot(211)
mesh(Lmatrix, Amatrix, v')
hold on
xlabel( ' labour productivity ', 'Fontsize', 14 )
ylabel( ' assets ','Fontsize', 14  )
title ( [' the value function with (\phi, \rho, \sigma) = (', num2str(phi), ', ', num2str(rhol), ', ', num2str(stdl),')  '],'Fontsize', 14  )
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

subplot(212)
mesh(Lmatrix, Amatrix, (g - phi*ones(size(g)))')
hold on
title ( ' decision rules ', 'Fontsize', 14  )
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

[Lmatrix, Amatrix2] = meshgrid(lgrid, muagrid - phi*ones(size(muagrid)));

figure

indicate = mu <= precision;
muse = mu; muse(indicate) = NaN;

mesh(Lmatrix, Amatrix2, muse')
hold on
xlabel( ' labour productivity ', 'Fontsize', 14 )
ylabel( ' assets ','Fontsize', 14  )
title ( [ ' stationary distribution when (\phi, \rho, \sigma) = (', num2str(phi), ', ', num2str(rhol), ', ', num2str(stdl),')  with  r = ' num2str(realrate) ' and K = ' num2str(kagg) ], 'Fontsize', 14  )
view(-150,20)
rotate3d on
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

figure

indicate = cons == 0;
consuse = cons;
consuse(indicate) = NaN;

mesh(Lmatrix, Amatrix2, consuse)
hold on
xlabel( ' labour productivity ', 'Fontsize', 14 )
ylabel( ' assets ','Fontsize', 14  )
title ( [ ' consumption with  r = ' num2str(realrate) ], 'Fontsize', 14  )
view(-150,20)
rotate3d on
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

figure

indicate = mpc == 0;
mpcuse = mpc;
mpcuse(indicate) = NaN;

mesh(Lmatrix, Amatrix2, mpcuse)
hold on
xlabel( ' labour productivity ', 'Fontsize', 14 )
ylabel( ' assets ','Fontsize', 14  )
title ( [ ' mpc when (\phi, \rho, \sigma) = (', num2str(phi), ', ', num2str(rhol), ', ', num2str(stdl),')  with  r = ' num2str(realrate) ' and K = ' num2str(kagg) ], 'Fontsize', 14  )
view(-150,20)
rotate3d on
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)


figure; 
yyaxis left; plot(muagrid, mpcall./muall, 'LineWidth', 3)
ylabel ( ' marginal propensity to consume ' )
yyaxis right; plot(muagrid, 100*muall, '.','LineWidth',3)
ylabel ( ' percentage of people ' )

title ([ ' aggregate mpc ' num2str(meanmpc) ] )
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

% new way to do elapsed time
finish = datetime('now');
elapsedtime = seconds(finish - starttime);

savingsrate = delta*alpha/(r + delta);

disp ([ ' elapsed time: ' num2str(elapsedtime) ])

% save results in a mat file that does not overwrite existing files and 
filecount = 0;
filexist = exist([resultfile '.mat'], 'file');
filename = resultfile;

while (filexist ~= 0)
     filecount = filecount + 1;
     filename = [resultfile '_' num2str(filecount) ];
     filexist = exist([filename '.mat'], 'file');
end

daytime = finish;

rss = r; wagess = wage; kvalss = kval; kaggss = kagg; laggss = lagg;
vss = v; gss = g; realratess = realrate; muss = mu; 

eval ( [ ' save ', filename, ...
          ' termsarray alpha delta pilr alow ahigh kl kh kvalss rss wagess kaggss laggss vss gss Lmatrix Amatrix Amatrix2 realratess muss savingsrate daytime' ] )

disp ( [ ' results saved as ', filename,'.mat ' ] ) 

disp ( ' ' )
disp ( ' ************************************************************* ')
disp ( ' norm in kagg is important for convergence of the distribution ' )
disp ( ' anum must be large enough, and agrid(anum) not too large, ' )
disp ( ' and precision small, to have kval - kagg small ' )

