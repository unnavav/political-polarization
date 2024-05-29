classdef egm
    methods(Static)

        %% getting derivative of expected value function
        % inputs: EV grid
        % outputs: DEV grid
        function DEV = numdev(EV, agrid)
            [nl, na, np] = size(EV);
            
            DEV = zeros(size(EV));

            for ip = 1:np
                for ia = 1:na
                    evia = EV(:, ia,ip);
                    a = agrid(ia);
    
                    if ia > 1 && ia < na
                        evia_0 = EV(:, ia-1,ip);
                        evia_1 = EV(:, ia+1,ip);
                        a1 = agrid(ia+1);
                        a0 = agrid(ia - 1);
                        
                        dr = (evia_1 - evia)/(a1 -a);
                        dl = (evia - evia_0)/(a -a0);
    
                        DEV(:, ia,ip) = (dr + dl)/2;
    
                    elseif ia == 1
    
                        evia_1 = EV(:, ia+1,ip);
                        a1 = agrid(ia+1);
                        dr = (evia_1 - evia)/(a1 -a);
    
                        DEV(:, ia,ip) = dr;
    
                    else
    
                        evia_0 = EV(:, ia-1,ip);
                        a0 = agrid(ia-1);
                        dl = (evia - evia_0)/(a - a0);
    
                        DEV(:, ia,ip) = dl;
    
                    end
    
                end
            end
        end

        function Dev = aubdev(ev,agrid)
            [m, anum] = size(ev);

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
        end

       function d = solveD(EV0, ik, kchgrid)
            nk = length(kchgrid);
            k = kchgrid(ik);

            if ik == nk
                d = (EV0(ik) - EV0(ik-1))/(k - kchgrid(ik-1));
            else
                d = (EV0(ik+1) - EV0(ik))/(kchgrid(ik+1) - k);
            end
        end

    end
end
