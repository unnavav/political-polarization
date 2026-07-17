function [zvec, piz] = tauchen(meanz, stdinnov, rho, multiple, znum)
   
% G. Tauchen (1986) 'Finite State Markov-Chain Approximations to 
% Univariate and Vector Autoregressions' Economics Letters 20: 177-181 
%
% [zvec, piz] = tauchen(meanz, stdinnov, rho, multiple, znum)
%
%               meanz       mean of stochastic process z
%               stdinnov    standard deviation of innovations to z
%               rho         persistence of stochastic process z
%               multiple    numerical parameter that determines bounds as multiples of the standard deviation of z
%               znum        the number of grid points in the discretised space
%               zvec        the discretised space
%               piz         Markov transition matrix on zvec, piz(i,j) = Pr{z' = z(j) | z = z(i)}

% bounds for tauchen grid are multiples of standard deviation of the
% stochastic process.
stdz = stdinnov^2.0;
stdz = stdz/(1.0 - rho^2.0);
stdz = sqrt(stdz);
zlow = meanz - stdz*multiple;
zhigh = meanz + stdz*multiple;
meaninnov = meanz*(1.0 - rho);


zvec = linspace(zlow, zhigh, znum);
piz = zeros(znum, znum);

w = (zhigh - zlow)/(znum-1);

% This is equations (3a) and (3b) from Tauchen (1986)

disp([zvec meaninnov stdinnov])

for j = 1:1:znum
    
    z = zvec(1) - rho*zvec(j);
    F1 = normal(z + w/2.0, meaninnov, stdinnov);
    piz(j,1) = F1;
%     fprintf("Row %i: F1 = %1.4f; mean = %2.4f; std = %2.4f\n", ...
%         j, F1, meaninnov, stdinnov)
    
    for k = 2:1:znum - 1
        z = zvec(k) - rho*zvec(j);
        F1 = normal(z + w/2.0, meaninnov, stdinnov);
        F0 = normal(z - w/2.0, meaninnov, stdinnov);
        piz(j,k) = F1 - F0;
    end
    
    z = zvec(znum) - rho*zvec(j);
    F0 = normal(z - w/2.0, meaninnov, stdinnov);
    piz(j,znum) = 1.0 - F0;
    
end


function [normal] = normal(z, mean, sd)

% Aludaat and Alodat (2008) /A Note on Approximating the Normal Distribution?
% Applied Mathematical Sciences 2(9), pages 425-29,

% convert to standard normal
x = (z - mean)/sd;

y = 0.7988*x*(1.0 + 0.04417*(x^2));

a = exp(2*y);
normal = a/(1.0 + a);
