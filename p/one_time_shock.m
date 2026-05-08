%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A ONE TIME CAPITAL SHOCK
% 2026 FEB
% VAASAVI
%
% DID YOU THINK I'D MAKE IT THIS FAR? BECAUSE i DIDN'T, I DIDN'T AT ALL.
% BUT READ IT AND WEEP. I'M NOW GOING TO SIMULATE A ONE-TIME CAPITAL
% DEPRECIATION SHOCK, SEND THIS TO MIDWEST MACRO, AND CRY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clear all; clc;
addpath(genpath(pwd));

%% pulling in steady states of interest and relevant policies

cd ../d/
load ks_rfore_endo_all.mat

cd ../p

Kpred = ks.forecastK(Kfore, Kgrid);
[EV1, EV2, Votes_EV] = gov.getVotingExpectations(V, pil, pid, Kpred, Kgrid);

%% generating data

T = 5000;
ddata = ones(T,1);
t0 = 2000;                 % start of crisis
len = 10;                  % 10 periods of high δ
ddata(t0:(t0+len-1)) = 2;  % force high depreciation


% get inital conditions

[pil, lgrid] = compute.getTauchen(nl,  mu, sigx, rho, range);

Kprdata = zeros(T,1);
distr_array = cell(T,1);
g0 = g1;
g0_cond = compute.condense(g0, amu, agrid);
K0 = dot(squeeze(sum(sum(g0_cond,1),2)),agrid);
Kprdata(1) = K0; distr_array{1} = g0;

% start with regime one, then check
Rdata = ones(T,1);
Prdata = Rdata;

for t = 2:1:T-1
    Kt = Kprdata(t-1);
    Rt = Rdata(t);
    dt = ddata(t);
    g_prev = distr_array{t-1};

    [ix, we] = compute.weight(Kgrid, Kt);
    g_t = we*G(ix, Rt, dt, :,:) + (1-we)*G(ix+1, Rt, dt, :, :);
    g_t = squeeze(g_t);
    g_today = HH.transitDistr(g_t, g_prev, amu, agrid, pil);

    distr_array{t} = g_today;
    acond = compute.condense(g_today, amu, agrid);
    Kpr = dot(squeeze(sum(sum(acond,1),2)),agrid);
    Kprdata(t) = Kpr;

    % now use K today, Kpr, and R today to back out max vote,
    % along with the actual distribution over wealth. Note that
    % this voting rule already interpolates over the future
    % forecast, so I only need to interpolate over today's K,
    % then force it back to binary (otherwise there's a decimal
    % value on whether or not I'll vote for R = 1)

    todays_votes = we*Votes_EV(ix, Rt, dt, :, :) + ...
        (1-we)*Votes_EV(ix+1, Rt, dt, :, :);
    todays_votes = squeeze(todays_votes);
    todays_votes = (todays_votes >= .5);

    vote_total = sum(sum(sum(squeeze(acond(dt,:,:)).*todays_votes)))*2;
    %multiplying by 2 bc each dimension has 50% of mass
    Prdata(t) = vote_total;
    if (vote_total <=.5) 
        Rdata(t+1) = 2; 
    else 
        Rdata(t+1) = 1;
    end

    if mod(t,100) == 0
        fprintf("\n\t t = %i", t)
    end
end

%% plotting

T = length(Kprdata);
t = (1900:2500)';

% Colors
pop_color  = [0.85 0.3 0.3];   % populism (R=1)
lib_color  = [0.3 0.3 0.7];    % liberalism (R=2)
cap_color  = [0 0 0];          % capital path

figure; hold on;

% 1) Shade regimes as background patches
ymin = min(Kprdata);
ymax = max(Kprdata);

% Find contiguous segments with same regime
R = Rdata(:);
start_idx = 1;
for tt = 2:T
    if R(tt) ~= R(tt-1)
        % draw block for [start_idx, tt-1]
        x0 = t(start_idx);
        x1 = t(tt-1);
        if R(start_idx) == 1
            c = pop_color;
        else
            c = lib_color;
        end
        patch([x0 x1 x1 x0],[ymin ymin ymax ymax],c, ...
            'FaceAlpha',0.12,'EdgeColor','none');
        start_idx = tt;
    end
end
% last segment
x0 = t(start_idx);
x1 = t(end);
if R(start_idx) == 1
    c = pop_color;
else
    c = lib_color;
end
patch([x0 x1 x1 x0],[ymin ymin ymax ymax],c, ...
    'FaceAlpha',0.12,'EdgeColor','none');

% 2) Plot capital on top
plot(t,Kprdata,'Color',cap_color,'LineWidth',1.6);

xlabel('Time');
ylabel('Aggregate capital K_t');
title('Capital Path with Populist (red) and Liberal (blue) Regimes');
set(gca,'FontSize',12);
box on;

%% try two
t    = (1:length(Kprdata))';   % full time index


start_idx = 1900;
end_idx   = 2500;

t_raw   = (start_idx:end_idx)';              % 1900:2500
t_shift = t_raw - 2000;                      % now runs from -100 to 500
Kall = Kprdata(:);
Rall = Rdata(:);

figure('Color','w'); hold on;

% y-range only over the window
Kwin = Kall(start_idx:end_idx);
ymin = min(Kwin);
ymax = max(Kwin);

pop_color  = [0.85 0.3 0.3];
lib_color  = [0.3 0.3 0.7];
cap_color  = [0 0 0];

% ---- shade regimes using FULL indices ----
reg_change = find(diff(Rwin) ~= 0);
seg_starts = [1; reg_change+1];
seg_ends   = [reg_change; length(Rwin)];

for s = 1:length(seg_starts)
    k0 = seg_starts(s);
    k1 = seg_ends(s);

    x0 = t_shift(k0);   % use shifted time
    x1 = t_shift(k1);

    if Rwin(k0) == 1
        c = pop_color;
    else
        c = lib_color;
    end

    patch([x0 x1 x1 x0],[ymin ymin ymax ymax],c, ...
          'FaceAlpha',0.12,'EdgeColor','none');
end

% Optional: purple band for 2042–2054
i0 = 2042 - start_idx + 1;   % index in window
i1 = 2054 - start_idx + 1;
x0p = t_shift(i0);
x1p = t_shift(i1);
patch([x0p x1p x1p x0p],[ymin ymin ymax ymax],[0.6 0.3 0.6], ...
      'FaceAlpha',0.18,'EdgeColor','none');
purple = [0.6 0.3 0.6];
x0 = 2042;
x1 = 2054;
patch([x0 x1 x1 x0],[ymin ymin ymax ymax],purple, ...
    'FaceAlpha',0.18,'EdgeColor','none');

% ---- capital path, restricted to window via xlim ----
plot(t_shift, Kwin, 'Color', cap_color, 'LineWidth', 1.6);
xlim([-100 500]);

xlabel('Time');
ylabel('Aggregate capital K_t');
title('Capital Path with Populist (red) and Liberal (blue) Regimes');
set(gca,'FontSize',12,'Color','w');
box on;
