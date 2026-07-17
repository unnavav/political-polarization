% aiyagimp transitional dynamics of Aiyagari model
%
%
% Boppart, Krusell and Mitman, 'Exploiting MIT shocks in
% heterogeneous-agent economies: the impulse response as a numerical 
% derivative' Journal of Economic Dynamics & Control 89 (2018) 68–92
%
% aiyagimpsim (this program) aiyagimprw (function computes perfect
% foresight path for kagg using inputs of rvec and wagevec)
%
% Aubhik Khan 2023

clear
close all
clc

% steady state file name
filename =  'resegmagphi1';

eval ( [ 'load ' filename ])

% filename contains the structure termsarray as well as alpha delta pilr alow ahigh kl kh kvalss rss wagess kaggss laggss vss gss Lmatrix Amatrix Amatrix2 realratess muss savingsrate daytime' ] )

% lnum = termsarray.lnum; lgrid = termsarray.lgrid; pil =  termsarray.pil; beta = termsarray.beta; sigma = termsarray.sigma; 
% anum = termsarray.anum; agrid = termsarray.agrid; phi = termsarray.phi; muagnum = termsarray.muagnum; muagrid = termsarray.muagrid;
% precision = termsarray.precision; precision1 = termsarray.precision1; 

% length of impulse and simulation
tmax = 150;
tsim = 1250;


% one time shock for impulse
innova0 = 0.01;

% persistence and std dev.
rhoa = 0.9;
stdinnova = 0.01;


% weight on existing guess for kvalin that sets r and w
lambda = 0.9;


atfpvec = ones(1, tmax);
atfpvec(1) = innova0;
for t = 1:1:tmax-1
    atfpvec(t+1) = rhoa*atfpvec(t);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterate on perfect foresight path for equilbrium capital stock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = exist('kvalimp.mat', 'file');

if (s == 2)
    load kvalimp
else
    kvalimp = ones(1, tmax-1)*kaggss;
end

% new way to time program
starttime = datetime('now');

lnow = laggss;
rvec = zeros(size(kvalimp));
wagevec = rvec;

% precision from structure
precision = termsarray.precision;

distance = 2.0*precision;
s1 = 0;

while (distance > precision/10 )


    for t = 1:1:tmax-1
        know = kvalimp(t);
        atfpval = exp(atfpvec(t));
        rvec(t) = alpha*atfpval*(know^(alpha - 1.0))*(lnow^(1.0 - alpha)) - delta;
        wagevec(t) = (1.0 - alpha)*atfpval*(know^alpha)*(lnow^(-alpha));
    end

    [kaggvec] = hw3q2p2rw(termsarray, tmax, rvec, wagevec, muss, vss);

    s1 = s1 + 1;
    distance = max(abs(kvalimp - kaggvec));

    fprintf( ' impulse %4d  innova0 = %8.4f   r(5)/rss = %8.2e   wage(5)/wss = %8.2e   distance = %8.2e  \n', s1, innova0, ...
                                                        (rvec(5) - rss)*100, 100*(wagevec(5) - wagess)/wagess, distance);
    disp ( ' ' )

    tkvalimp = lambda*kvalimp + (1.0 - lambda)*kaggvec;

    kvalimp = tkvalimp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impulse response as derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dkvalimp = (kaggvec - kaggss)/innova0;

figure('units','normalized','outerposition',[0 0 1 1])
plot(dkvalimp, 'LineWidth',3); 
legend(' dkvalimp ')
title ( ' capital response to innovation $\frac{dk_{t}}{da_{0}}$ ', 'Interpreter', 'latex')
set(gca,'FontSize',20);

% dkvalimp(s) is the derivative of k(t+s) for a shock to tfp at t.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save kvalimp kvalimp

rng(886644, 'twister')

innovec = stdinnova*randn(1,tsim);
% innovec = zeros(1, tsim);
% innovec(1) = 0.01;

atfplag = 0; atfpsim = zeros(1, tsim);
for t = 1:1:tsim
    atfpsim(t) = rhoa*atfplag + innovec(t);
    atfplag = atfpsim(t);
end

atfpsim = exp(atfpsim);
dksim = zeros(1, tsim);


% dk(1) = 0, dkvalimp(1) = 0, capital is a state
% dk(2) = dkvalimp(2)*innov(1)
% dk(3) = dkvalimp(3)*inov(1) + dkvalimp(2)*innov(2)
% dk(4) = dkvalimp(4)*innov(1) + dkvalimp(3)*innov(2) + dkvalimp(2)*innov(3)

% Take tmax = 5, then dk(1) to dk(4) are given above, 
% dk(5) = dkvalimp(4)*innov(2) + dkvalimp(3)*innov(3) + dkvalimp(2)*innov(4)

% impulse has dkvalimp(1),...,dkvalimp(tmax-1)
% dksim(t+1), for t = 1,...,tmax-1 will depend on innov(1),...,innov(t)
% dksim(t+1) = dkvalimp(2)*innov(t) + dkvalimp(3)*innov(t-1) + ...
% dkvalimp(t)*innov(1).  
% when t = tmax,...tsim
% dksim(t+1) = dkvalimp(2)*innov(t) + dkvalimp(3)*innov(t-1) + ...
% dkvalimp(tmax-1)*innov(t - (tmax-2))

% version from leisure model, debugged at 1559 on 12 November 23.
dksim2 = dksim;

for t = 1:1:tsim - 1
    for s = 1:min(t-1, tmax-2)
        dksim2(t+1) = dksim2(t+1) + dkvalimp(s+1)*innovec(t + 1 - s);
    end
end
ksim2 = dksim2 + kaggss;

ksimu = ksim2;

ysim = atfpsim.*(ksimu.^alpha).*(lnow.^(1.0-alpha));
isim = ksimu(2:tsim) - (1.0 - delta)*ksimu(1:tsim-1);
isim(tsim) = delta*kaggss;
csim = ysim - isim;

iaggss = delta*kaggss;
yaggss = (kaggss.^alpha).*(lnow.^(1.0-alpha));
caggss = yaggss - iaggss;

chat = (csim - caggss)/caggss; chat = 100*chat;
ihat = (isim - iaggss)/iaggss; ihat = 100*ihat;
yhat = (ysim - yaggss)/yaggss; yhat = 100*yhat;

% new way to do elapsed time
finish = datetime('now');
elapsedtime = seconds(finish - starttime);

disp ([ ' elapsed time: ' num2str(elapsedtime) ])


time = tmax:1:tsim;

figure('units','normalized','outerposition',[0 0 1 1])
subplot(211)
plot(time, atfpsim(time), 'LineWidth',3); 
legend(' tfp ', 'Location', 'southeast')
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

subplot(212)
plot(time, ksimu(time), 'LineWidth',3); 
legend (' k ', 'Location', 'southeast' )
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

figure('units','normalized','outerposition',[0 0 1 1])
subplot(311)
plot(time, yhat(time), 'LineWidth',3); 
legend (' y ', 'Location','southeast' )
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

subplot(312)
plot(time, chat(time), 'm','LineWidth',3); 
legend (' c ', 'Location','southeast' )
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

subplot(313)
plot( time, ihat(time), 'LineWidth',3); 
legend (' inv ', 'Location','southeast' )
set(gca,'FontSize',20);
set(get(gca,'Xlabel'),'FontSize',20)

