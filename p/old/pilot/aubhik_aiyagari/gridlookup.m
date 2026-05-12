function [iloc] = gridlookup(n, grid, valplace)

% Aubhik, Dublin, 31 1 2009.  

% This function will return iloc such that grid(iloc) is the greatest 
% value less than or equal to valplace, unless valplace is less than 
% grid(1) in which case iloc = 1.  

ilow = 1; ihigh = n;

distance = 2;

while (distance > 1)
    
    inow = floor((ilow + ihigh)/2);
    valnow = grid(inow);
    
    % The strict inequality here ensures that grid(iloc) is less than or
    % equal to valplace.
    if (valnow > valplace)
        ihigh = inow;
    else
        ilow = inow;
    end
    
    distance = ihigh - ilow;
    
end

iloc = ilow;