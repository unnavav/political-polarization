classdef egm
    methods(Static)

        %% getting derivative of expected value function
        % inputs: EV grid
        % outputs: DEV grid
        function DEV = numdev(EV, agrid)
            [nl, na, np] = size(EV);
            
            DEV = zeros(size(EV));

            for ip = 1:np
                for il = 1:nl
                    for ia = 1:na
        
                        evia = EV(il, ia,ip);
                        a = agrid(ia);
        
                        if ia > 1 && ia < na
                            evia_0 = EV(il, ia-1,ip);
                            evia_1 = EV(il, ia+1,ip);
                            a1 = agrid(ia+1);
                            a0 = agrid(ia - 1);
                            
                            dr = (evia_1 - evia)/(a1 -a);
                            dl = (evia - evia_0)/(a -a0);
        
                            DEV(il, ia,ip) = (dr + dl)/2;
        
                        elseif ia == 1
        
                            evia_1 = EV(il, ia+1,ip);
                            a1 = agrid(ia+1);
                            dr = (evia_1 - evia)/(a1 -a);
        
                            DEV(il, ia,ip) = dr;
        
                        else
        
                            evia_0 = EV(il, ia-1,ip);
                            a0 = agrid(ia-1);
                            dl = (evia - evia_0)/(a - a0);
        
                            DEV(il, ia,ip) = dl;
        
                        end
        
                    end
                end
            end
        end

    end
end
