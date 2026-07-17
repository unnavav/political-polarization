% perfect foresight transitional dynamics for the one sector growth model

clear
clc
close all

% parameters are capital's share, the depreciation rate for capital, the
% subjective discount factor and the inverse of the elasticity of
% intertemporal substitution of consumption.  
alpha = 0.36;
delta = 0.06;
beta = 0.96;
sigma = 2.0;

% Additional parameters is the length of the transition.  
simlength = 175;

% smoothing parameters which governs the speed with which the path of
% capital is updated.  
lambda = 0.6; 

% steady state of the model
zss = 1.0;
kterm = 1.0/beta - (1.0 - delta);
kss = alpha*zss/kterm;
kss = kss^(1.0/(1.0 - alpha));
css = zss*(kss^alpha) - delta*kss;

yss = zss*(kss^alpha);
iss = delta*kss;
rss = alpha*yss/kss - delta;

disp ( ' Perfect Foresight Transitional Dynamics of the One-Sector Growth Model ' )
disp ( ' ' )
disp ( ' remark: if T is too small, then the model will not return to the steady state ' )
disp ( ' ' )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%      Transitional Dyanmics        %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


zvec = ones(1,simlength-1)*zss;
kvec = ones(1, simlength-1)*kss;

precision = 0.0001;
distance = 2.0*precision;
iteration = 0;

kvec(1) = 0.1*kss;

while (distance > precision)
    
    % solve for real interest rates implied by the conjectured path of k
    rvec = alpha*zvec.*(kvec.^(alpha - 1.0)) - delta;
      
    % solve for path of consumption, critically note that if r(T+1) is the
    % steady state level, then c(T-1) must equal css.  Starting at this
    % value, we solve backwards for the path of consumption implied by real
    % interest rates.
    cvec = zeros(1, simlength - 1);
    
    t = simlength - 1;
    cvec(t) = css;
    for t = simlength - 2:-1:1
        eulerterm = beta*(1.0 + rvec(t+1));
        eulerterm = eulerterm^(-1.0/sigma);
        cvec(t) = eulerterm*cvec(t+1);
    end
    
    % Output in each period <i>using the conjectured path for k<i> is used,
    % alongside consumption, to derive an implied path that is Tk.
    yvec = zvec.*(kvec.^alpha);
    
    Tkvec = zeros(1, simlength - 1);
    Tkvec(1) = kvec(1);
    for t = 1:1:simlength - 2
        Tkvec(t+1) = (1.0 - delta)*Tkvec(t) + yvec(t) - cvec(t);
    end
        
    % The algorithm stops when k = Tk.  
    distance = max(abs(kvec - Tkvec)); iteration = iteration + 1;
    s = sprintf ( ' iteration = %4d     ||Tk - k|| = %8.4f ', iteration, distance);
    disp(s)
    
    kvec = lambda.*kvec + (1.0 - lambda).*Tkvec;

end

ivec = [kvec(2:simlength-1) kss] - (1.0 - delta)*kvec(1:simlength-1);

time = 1:1:floor((simlength-1)/3);
figure

subplot(311)
plot(time, kvec(time), time, ones(size(time))*kss, 'LineWidth',4)
legend( ' k_{t} ', ' k_{steady state} ', 'Location', 'SE')
ylabel ( ' level ', 'FontSize', 14 )

subplot(312)
plot(time, yvec(time), time, ones(size(time))*yss, 'LineWidth',4)
legend(' y_{t} ', ' y_{steady state} ', 'Location', 'SE')
ylabel ( ' level', 'FontSize', 14 )

subplot(313)
plot(time, cvec(time), time, ones(size(time))*css, 'LineWidth',4)
legend( ' c_{t} ', ' c_{steady state} ', 'Location', 'SE')
ylabel ( ' level ', 'FontSize', 14 )

    
    

