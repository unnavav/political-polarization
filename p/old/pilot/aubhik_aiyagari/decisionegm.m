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