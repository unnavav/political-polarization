classdef egm
    methods(Static)

        %% getting derivative of expected value function
        % inputs: EV grid
        % outputs: DEV grid
        function DEV = numdev(EV, agrid)
            [nl, na] = size(EV);
            
            DEV = zeros(size(EV));

            for il = 1:nl
                for ia = 1:na
                    evia = EV(il, ia);
                    a = agrid(ia);
    
                    if ia > 1 && ia < na
                        evia_0 = EV(il, ia-1);
                        evia_1 = EV(il, ia+1);
                        a1 = agrid(ia+1);
                        a0 = agrid(ia - 1);
                        
                        dr = (evia_1 - evia)/(a1 -a);
                        dl = (evia - evia_0)/(a -a0);
    
                        DEV(il, ia) = (dr + dl)/2;
    
                    elseif ia == 1
    
                        evia_1 = EV(il, ia+1);
                        a1 = agrid(ia+1);
                        dr = (evia_1 - evia)/(a1 -a);
    
                        DEV(il, ia) = dr;
    
                    else
    
                        evia_0 = EV(il, ia-1);
                        a0 = agrid(ia-1);
                        dl = (evia - evia_0)/(a - a0);
    
                        DEV(il, ia) = dl;
    
                    end
    
                end
            end
        end

        %% gives us the function of the budget constraint so we can 
        % find the implied a
        % inputs: c+apr, w*eps, r, tau
        % outputs: value of budget constraint for given x
        function budget = impliedbc(x, bcterms)
            capr = bcterms.capr;
            we = bcterms.we;
            r = bcterms.r;
            tau = bcterms.tau;

            budget = (1+r)*x + we - (r*x + we)^(1-tau) - capr;

        end
    end
end
