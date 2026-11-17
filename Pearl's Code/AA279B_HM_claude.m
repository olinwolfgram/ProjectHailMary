%% AA279B_HM_Merged.m
% Earth -> Venus (powered flyby) -> Earth trajectory search, 2020-2030
% Patched conics, sun-only heliocentric legs, Curtis Lambert (prograde, 0 rev)
%
% Stage 0   Constants, settings, ephemeris precompute
% Stage 1a  Coarse scan of Earth -> Venus launch C3 to locate launch windows
% Stage 1   Fine Earth -> Venus Lambert grid per window (launch C3, DLA, Venus V-inf in)
% Stage 2   Fine Venus -> Earth Lambert grid per window (Venus V-inf out, Earth arrival C3)
% Stage 3   Venus flyby: feasibility filter, lower-bound pruning, r_p via fzero,
%           dV_flyby, cost J (Option B)
% Stage 4   Selection: optimal trajectory per launch period (2020-2023, 2023-2026,
%           2026-2030), plus the best <= 1 yr trajectory in each period if it differs
% Stage 5   ode113 propagation of the reported trajectories (sun gravity only)
% Stage 6   One plot per trajectory with labeled maneuvers, delta-V and dates;
%           summary table
%
% External dependencies (unchanged): AA279lambert_curtis,
% AA279j2000_planetary_elements, rot

clear; close all; clc;

%% Stage 0: Constants and settings ------------------------------------------

% Grav parameter (km^3/s^2)
mu.venus = 324858.630;
mu.earth = 398600.440;
mu_sun = 1.3271244004193938e11;

% Mean radius (km)
r_m.venus = 6051.80;
r_m.earth = 6371.01;

AU = 149597870.7;              % km
eps_J2000 = deg2rad(23.43928); % obliquity at J2000; ephemeris assumed J2000 ecliptic

% Venus flyby limits
h_min = 300;                   % km, minimum flyby altitude (mission choice)
r_min = r_m.venus + h_min;     % km
[el_venus, ~] = AA279j2000_planetary_elements(2);
r_SOI = el_venus(1)*AU * (mu.venus/mu_sun)^(2/5); % km, Venus SOI (~6.2e5 km)

% Arrival and Departures
p1_ID = 3; % Earth
p2_ID = 2; % Venus

t_begin = mjuliandate(datetime([2020, 01, 01, 00, 00, 00], 'TimeZone','UTC'));
t_end   = mjuliandate(datetime([2030, 12, 31, 00, 00, 00], 'TimeZone','UTC'));

step = 2;                          % days, fine-grid step
step_coarse = 10;                  % days, window-finding step
arr_times = 100 : step : 150;      % days, Earth -> Venus TOF
return_times = 80 : step*2 : 350;  % days, Venus -> Earth TOF

% Launch window finding
win_sep = 400;        % days, min spacing between windows (synodic period ~584 d)
win_half = 45;        % days, fine-grid half width around each window center
C3_window_max = 30;   % km^2/s^2, ignore coarse minima above this launch C3

% Selection: launch periods [2020, 2023), [2023, 2026), [2026, end of 2030]
t_2023 = mjuliandate(datetime([2023, 01, 01, 00, 00, 00], 'TimeZone','UTC'));
t_2026 = mjuliandate(datetime([2026, 01, 01, 00, 00, 00], 'TimeZone','UTC'));
period_edges = [t_begin, t_2023, t_2026, t_end];
period_names = ["2020-2023", "2023-2026", "2026-2030"];
n_period = numel(period_names);
T_max = 365;          % days, total duration limit for the "best <= 1 yr" option

options = optimset('TolX', 1e-3, 'Display', 'off'); % km

mjd2str = @(m) string(datetime(m, 'ConvertFrom', 'modifiedjuliandate', ...
                               'Format', 'dd-MMM-yyyy'));

%% Stage 1a: Coarse scan to locate launch windows -----------------------------
dep_dates_c = t_begin : step_coarse : t_end;
arr_times_c = arr_times(1) : step_coarse : arr_times(end);

[X_E_c, V_E_c] = eph_states(p1_ID, dep_dates_c);
arr_grid_c = dep_dates_c(:) + arr_times_c;
[arr_dates_c, ~, ic] = unique(arr_grid_c(:));
[X_V_c, V_V_c] = eph_states(p2_ID, arr_dates_c);
arr_idx_c = reshape(ic, size(arr_grid_c));

p1_Vinf_c = lambert_grid(X_E_c, V_E_c, X_V_c, V_V_c, arr_idx_c, arr_times_c, mu_sun);
p1_C3_c = reshape(sum(p1_Vinf_c.^2, 1), numel(dep_dates_c), numel(arr_times_c));

C3_c_min = min(p1_C3_c, [], 2);   % best launch C3 per departure date
is_win = islocalmin(C3_c_min, 'SamplePoints', dep_dates_c(:), ...
                    'MinSeparation', win_sep) & (C3_c_min <= C3_window_max);
win_centers = dep_dates_c(is_win);
if isempty(win_centers)
    error('No launch windows found; check C3_window_max or the scan range.');
end
fprintf('Launch windows found: %d\n', numel(win_centers));

%% Stages 1-3: fine grid in each window ---------------------------------------
cand_mat = zeros(0, 13);

% Running best J per launch period (all durations, and <= T_max), used to prune
J_best_all = Inf(1, n_period);
J_best_1yr = Inf(1, n_period);

for w = 1:numel(win_centers)
    dep_dates = max(t_begin, win_centers(w) - win_half) : step : ...
                min(t_end,   win_centers(w) + win_half);
    Nd = numel(dep_dates); Na = numel(arr_times); Nr = numel(return_times);

    % Stage 0: ephemerides only on the unique dates this window needs
    fb_grid = dep_dates(:) + arr_times;             % Nd x Na, Venus flyby MJD
    [fb_dates, ~, ic] = unique(fb_grid(:));
    fb_idx = reshape(ic, Nd, Na);
    Nfb = numel(fb_dates);
    ret_grid = fb_dates + return_times;             % Nfb x Nr, Earth return MJD
    [ret_dates, ~, ir] = unique(ret_grid(:));
    ret_idx = reshape(ir, Nfb, Nr);

    [X_E_dep, V_E_dep] = eph_states(p1_ID, dep_dates);
    [X_V_fb,  V_V_fb ] = eph_states(p2_ID, fb_dates);
    [X_E_ret, V_E_ret] = eph_states(p1_ID, ret_dates);

    % Stage 1: Earth -> Venus
    [p1_Vinf_grid, p2_Vinf_grid] = lambert_grid(X_E_dep, V_E_dep, X_V_fb, V_V_fb, ...
                                                fb_idx, arr_times, mu_sun);
    p1_C3_total_grid = reshape(sum(p1_Vinf_grid.^2, 1), Nd, Na);  % launch C3
    p2_C3_total_grid = reshape(sum(p2_Vinf_grid.^2, 1), Nd, Na);  % Venus arrival V-inf^2

    % Launch asymptote declination, rotated ecliptic -> Earth equatorial
    vz_eq = reshape(p1_Vinf_grid(2,:,:)*sin(eps_J2000) + ...
                    p1_Vinf_grid(3,:,:)*cos(eps_J2000), Nd, Na);
    DLA1_total_grid = asind(vz_eq ./ sqrt(p1_C3_total_grid));    % deg

    % Stage 2: Venus -> Earth
    [p1_Vinf_r_grid, p2_Vinf_r_grid] = lambert_grid(X_V_fb, V_V_fb, X_E_ret, V_E_ret, ...
                                                    ret_idx, return_times, mu_sun);
    p1_C3_return = reshape(sum(p2_Vinf_r_grid.^2, 1), Nfb, Nr);  % Earth arrival V-inf^2

    % Stage 3: Venus flyby. Rows = leg-1 cells (j fastest), cols = return TOF
    n1 = Nd*Na;
    vin = reshape(p2_Vinf_grid, 3, n1);             % incoming V-inf at Venus
    vix = vin(1,:).'; viy = vin(2,:).'; viz = vin(3,:).';
    Vox = reshape(p1_Vinf_r_grid(1,:,:), Nfb, Nr); Vox = Vox(fb_idx(:), :);
    Voy = reshape(p1_Vinf_r_grid(2,:,:), Nfb, Nr); Voy = Voy(fb_idx(:), :);
    Voz = reshape(p1_Vinf_r_grid(3,:,:), Nfb, Nr); Voz = Voz(fb_idx(:), :);

    % Required turn angle between incoming and outgoing V-inf (atan2 form)
    dotp = vix.*Vox + viy.*Voy + viz.*Voz;
    crs  = sqrt((viy.*Voz - viz.*Voy).^2 + (viz.*Vox - vix.*Voz).^2 + ...
                (vix.*Voy - viy.*Vox).^2);
    del_tot = atan2(crs, dotp);

    vin2  = repmat(p2_C3_total_grid(:), 1, Nr);
    vout2 = Vox.^2 + Voy.^2 + Voz.^2;

    % Feasibility filter: turn available at r_min (max) and r_SOI (min)
    delta_max = asin(1./(1 + r_min*vin2/mu.venus)) + asin(1./(1 + r_min*vout2/mu.venus));
    delta_min = asin(1./(1 + r_SOI*vin2/mu.venus)) + asin(1./(1 + r_SOI*vout2/mu.venus));

    dep_MJD = repmat(dep_dates(:), Na, Nr);
    fb_MJD  = repmat(fb_dates(fb_idx(:)), 1, Nr);
    ret_MJD = fb_dates(fb_idx(:)) + return_times;

    feasible = (del_tot <= delta_max) & (del_tot >= delta_min) & (ret_MJD <= t_end);

    % Exact lower bound on J without solving for r_p. dV_flyby increases with
    % r_p (its denominator sqrt(vout2+2mu/r)+sqrt(vin2+2mu/r) shrinks as r grows),
    % so dV_flyby(r_p) >= dV_flyby(r_min) for every feasible cell.
    C3_launch = repmat(p1_C3_total_grid(:), 1, Nr);
    C3_return = p1_C3_return(fb_idx(:), :);
    dV_lb = abs(sqrt(vout2 + 2*mu.venus/r_min) - sqrt(vin2 + 2*mu.venus/r_min));
    J_lb  = sqrt(C3_launch) + dV_lb + sqrt(C3_return);

    period = discretize(dep_MJD, period_edges);    % launch period 1..3
    is_1yr = (ret_MJD - dep_MJD) <= T_max;

    % Solve r_p only where the cell could still beat a running best (cheapest first)
    cells = find(feasible & isfinite(J_lb));
    [~, ord] = sort(J_lb(cells));
    cells = cells(ord);

    r_p = NaN(n1, Nr); dV_flyby = NaN(n1, Nr); J = NaN(n1, Nr);
    for idx = cells.'
        p = period(idx);
        if J_lb(idx) >= J_best_all(p) && ~(is_1yr(idx) && J_lb(idx) < J_best_1yr(p))
            continue
        end

        % Periapsis radius (bracket sign change guaranteed by the filter)
        r_p(idx) = fzero(@(x) fun(x, mu.venus, sqrt(vin2(idx)), sqrt(vout2(idx)), ...
                                  del_tot(idx)), [r_min r_SOI], options);

        % Periapsis delta-V of the powered flyby
        dV_flyby(idx) = abs(sqrt(vout2(idx) + 2*mu.venus/r_p(idx)) - ...
                            sqrt(vin2(idx)  + 2*mu.venus/r_p(idx)));

        % Cost (Option B), km/s
        J(idx) = sqrt(C3_launch(idx)) + dV_flyby(idx) + sqrt(C3_return(idx));

        J_best_all(p) = min(J_best_all(p), J(idx));
        if is_1yr(idx)
            J_best_1yr(p) = min(J_best_1yr(p), J(idx));
        end
    end

    DLA = repmat(DLA1_total_grid(:), 1, Nr);
    keep = find(isfinite(J));
    cand_mat = [cand_mat; J(keep), dep_MJD(keep), fb_MJD(keep), ret_MJD(keep), ...
                C3_launch(keep), DLA(keep), r_p(keep) - r_m.venus, dV_flyby(keep), ...
                sqrt(vin2(keep)), sqrt(vout2(keep)), sqrt(C3_return(keep)), ...
                w*ones(numel(keep), 1), period(keep)]; %#ok<AGROW>

    fprintf('Window %d (launch ~%s): %d feasible triplets, %d flyby solves\n', w, ...
            mjd2str(win_centers(w)), numel(cells), numel(keep));
end

if isempty(cand_mat)
    error('No feasible Earth-Venus-Earth trajectories found.');
end
cand = array2table(cand_mat, 'VariableNames', {'J','dep_MJD','fb_MJD','ret_MJD', ...
    'C3_launch','DLA','h_p','dV_flyby','vinf_fb_in','vinf_fb_out','vinf_arr','window', ...
    'period'});

%% Stage 4: Selection ----------------------------------------------------------
T_total = ret_MJD - cand.dep_MJD;                     % days

idx_report = zeros(1, 0);
tags = strings(1, 0);
for p = 1:n_period
    J_p = cand.J;  J_p(cand.period ~= p) = NaN;
    [J_opt, i_opt] = min(J_p);
    if isnan(J_opt)
        fprintf('No feasible trajectory with launch in %s.\n', period_names(p));
        continue
    end

    J_p1 = J_p;  J_p1(T_total > T_max) = NaN;
    [J_1, i_1] = min(J_p1);

    if isnan(J_1)
        idx_report(end+1) = i_opt;                                         %#ok<AGROW>
        tags(end+1) = period_names(p) + " optimal";                        %#ok<AGROW>
        fprintf('No feasible trajectory <= %d days with launch in %s.\n', ...
                T_max, period_names(p));
    elseif i_1 == i_opt
        idx_report(end+1) = i_opt;                                         %#ok<AGROW>
        tags(end+1) = period_names(p) + " optimal (also best <= 1 yr)";    %#ok<AGROW>
    else
        idx_report(end+1:end+2) = [i_opt, i_1];                            %#ok<AGROW>
        tags(end+1:end+2) = [period_names(p) + " optimal", ...
                             period_names(p) + " best <= 1 yr"];          %#ok<AGROW>
    end
end
if isempty(idx_report)
    error('No trajectories selected.');
end
n_traj = numel(idx_report);

%% Stage 5: Propagation with ode113 (sun gravity only) ------------------------
opts_ode = odeset('RelTol', 1e-12, 'AbsTol', [1e-3*ones(3,1); 1e-9*ones(3,1)]);
f_2b = @(t, y) [y(4:6); -mu_sun*y(1:3)/norm(y(1:3))^3];

traj = struct([]);
for q = 1:n_traj
    c = cand(idx_report(q), :);
    TOF1 = c.fb_MJD - c.dep_MJD;   % days
    TOF2 = c.ret_MJD - c.fb_MJD;   % days

    [x_E0, v_E0] = heliocentric(p1_ID, c.dep_MJD);
    [x_V,  v_V ] = heliocentric(p2_ID, c.fb_MJD);
    [x_E1, v_E1] = heliocentric(p1_ID, c.ret_MJD);

    [~, ~, v1_out,   v2_out] = p1_To_p2_dV(x_E0, v_E0, x_V, v_V, TOF1, mu_sun);
    [~, ~, v1_out_r, ~     ] = p1_To_p2_dV(x_V, v_V, x_E1, v_E1, TOF2, mu_sun);

    [~, Y1] = ode113(f_2b, linspace(0, TOF1*86400, 500), [x_E0; v1_out],   opts_ode);
    [~, Y2] = ode113(f_2b, linspace(0, TOF2*86400, 500), [x_V;  v1_out_r], opts_ode);

    traj(q).Y1 = Y1;
    traj(q).Y2 = Y2;
    traj(q).err1 = norm(Y1(end,1:3).' - x_V);    % km, leg-1 end vs Venus
    traj(q).err2 = norm(Y2(end,1:3).' - x_E1);   % km, leg-2 end vs Earth
    traj(q).dv_grav = norm(v1_out_r - v2_out);   % km/s, heliocentric change from flyby
end

%% Stage 6: Plots (one per trajectory) and summary ------------------------------
for q = 1:n_traj
    c = cand(idx_report(q), :);
    Y1 = traj(q).Y1(:,1:3)/AU;
    Y2 = traj(q).Y2(:,1:3)/AU;

    figure; hold on; grid on; axis equal;
    plot3(0, 0, 0, 'o', 'Color', [0.9 0.7 0], 'MarkerFaceColor', [1 0.85 0], ...
          'MarkerSize', 10, 'DisplayName', 'Sun');

    % Planet orbits sampled from this trajectory's launch date
    X_E_orb = eph_states(p1_ID, c.dep_MJD + linspace(0, 365.25, 361));
    X_V_orb = eph_states(p2_ID, c.dep_MJD + linspace(0, 224.70, 361));
    plot3(X_E_orb(1,:)/AU, X_E_orb(2,:)/AU, X_E_orb(3,:)/AU, ':', 'Color', [0 0.3 0.8], ...
          'DisplayName', 'Earth orbit');
    plot3(X_V_orb(1,:)/AU, X_V_orb(2,:)/AU, X_V_orb(3,:)/AU, ':', 'Color', [0.85 0.5 0], ...
          'DisplayName', 'Venus orbit');

    % Legs
    plot3(Y1(:,1), Y1(:,2), Y1(:,3), '-',  'Color', [0 0.45 0.74], 'LineWidth', 1.5, ...
          'DisplayName', 'Leg 1: Earth \rightarrow Venus');
    plot3(Y2(:,1), Y2(:,2), Y2(:,3), '-',  'Color', [0.64 0.08 0.18], 'LineWidth', 1.5, ...
          'DisplayName', 'Leg 2: Venus \rightarrow Earth');

    % Launch
    pL = Y1(1,:);
    plot3(pL(1), pL(2), pL(3), 'o', 'Color', 'k', 'MarkerFaceColor', [0 0.45 0.74], ...
          'MarkerSize', 8, 'DisplayName', 'Launch');
    text(pL(1), pL(2), pL(3), ...
        {sprintf(' Launch: %s', mjd2str(c.dep_MJD)), ...
         sprintf(' C3 = %.2f km^2/s^2, V_\\infty = %.2f km/s', c.C3_launch, sqrt(c.C3_launch)), ...
         sprintf(' DLA = %.1f deg', c.DLA)}, 'FontSize', 8);

    % Venus flyby
    pF = Y2(1,:);
    plot3(pF(1), pF(2), pF(3), 's', 'Color', 'k', 'MarkerFaceColor', [0.85 0.5 0], ...
          'MarkerSize', 8, 'DisplayName', 'Venus flyby');
    text(pF(1), pF(2), pF(3), ...
        {sprintf(' Venus flyby: %s', mjd2str(c.fb_MJD)), ...
         sprintf(' h_p = %.0f km, \\DeltaV_{flyby} = %.3f km/s', c.h_p, c.dV_flyby), ...
         sprintf(' |\\Deltav_{grav}| = %.2f km/s (gravity, not propulsive)', traj(q).dv_grav)}, ...
         'FontSize', 8);

    % Earth arrival
    pA = Y2(end,:);
    plot3(pA(1), pA(2), pA(3), '^', 'Color', 'k', 'MarkerFaceColor', [0.64 0.08 0.18], ...
          'MarkerSize', 8, 'DisplayName', 'Earth arrival');
    text(pA(1), pA(2), pA(3), ...
        {sprintf(' Earth arrival: %s', mjd2str(c.ret_MJD)), ...
         sprintf(' V_\\infty = %.2f km/s', c.vinf_arr), ...
         sprintf(' Total = %d days, J = %.3f km/s', round(T_total(idx_report(q))), c.J)}, ...
         'FontSize', 8);

    xlabel('x [AU]'); ylabel('y [AU]'); zlabel('z [AU]');
    title({sprintf('T%d: Earth \\rightarrow Venus flyby \\rightarrow Earth', q), ...
           char(tags(q))});
    legend('Location', 'bestoutside');
    view(2);
end

ir = idx_report(:);
Summary = table((1:n_traj).', tags(:), mjd2str(cand.dep_MJD(ir)), ...
    mjd2str(cand.fb_MJD(ir)), mjd2str(cand.ret_MJD(ir)), ...
    cand.fb_MJD(ir) - cand.dep_MJD(ir), cand.ret_MJD(ir) - cand.fb_MJD(ir), T_total(ir), ...
    cand.C3_launch(ir), sqrt(cand.C3_launch(ir)), cand.DLA(ir), cand.h_p(ir), ...
    cand.dV_flyby(ir), [traj.dv_grav].', cand.vinf_arr(ir), cand.J(ir), ...
    [traj.err1].', [traj.err2].', ...
    'VariableNames', {'Option','Tag','Launch','Venus_flyby','Earth_arrival', ...
    'TOF1_d','TOF2_d','Total_d','C3_launch_km2s2','Vinf_launch_kms','DLA_deg', ...
    'h_p_km','dV_flyby_kms','dv_grav_kms','Vinf_arr_kms','J_kms', ...
    'EndErr_Venus_km','EndErr_Earth_km'});
disp(Summary);

%% Functions -------------------
% Interplanetary Lambert wrapper (precomputed heliocentric states, dt in days)
function [p1_Vinf, p2_Vnif, v1_out, v2_out] = p1_To_p2_dV(x_dep_hc, v_dep_hc, ...
                                                          x_arr_hc, v_arr_hc, dt, mu_sun)
dm = 'pro';
nrev = 0;
[v1_out, v2_out, error_out] = AA279lambert_curtis(mu_sun, ...
                            x_dep_hc, x_arr_hc, dm, nrev, dt*24*3600);
v1_out = v1_out(:); v2_out = v2_out(:);

% AA279lambert_curtis reports "finished" on a successful solve
if numel(v1_out) ~= 3 || numel(v2_out) ~= 3 || error_out ~= "finished" || ...
        ~all(isfinite([v1_out; v2_out]))
    p1_Vinf = NaN(3,1); p2_Vnif = NaN(3,1);
    v1_out  = NaN(3,1); v2_out  = NaN(3,1);
    return
end

% Find planetocentric Vinf
p1_Vinf = v1_out - v_dep_hc;
p2_Vnif = v2_out - v_arr_hc;
end

% Lambert solves over (departure date j) x (TOF k)
% X1, V1: 3 x N1 departure-planet states; X2, V2: arrival-planet states at the
% unique arrival dates; idx2(j,k): column of X2/V2 for arrival at dep(j) + TOF(k)
function [Vinf_dep, Vinf_arr] = lambert_grid(X1, V1, X2, V2, idx2, TOF, mu_sun)
N1 = size(X1, 2); Nt = numel(TOF);
Vinf_dep = NaN(3, N1, Nt);
Vinf_arr = NaN(3, N1, Nt);
for j = 1:N1   % candidate for parfor (Parallel Computing Toolbox)
    for k = 1:Nt
        q = idx2(j, k);
        [Vinf_dep(:,j,k), Vinf_arr(:,j,k)] = p1_To_p2_dV(X1(:,j), V1(:,j), ...
                                             X2(:,q), V2(:,q), TOF(k), mu_sun);
    end
end
end

% Heliocentric states of one planet on a list of MJDs
function [X, V] = eph_states(p_ID, MJD_list)
n = numel(MJD_list);
X = zeros(3, n); V = zeros(3, n);
for q = 1:n
    [X(:,q), V(:,q)] = heliocentric(p_ID, MJD_list(q));
end
end

function [x_hc, v_hc] = heliocentric(p_ID, MJD)
    % Element tables loaded once per planet
    persistent elements_cache rates_cache
    if isempty(elements_cache)
        elements_cache = cell(9,1); rates_cache = cell(9,1);
    end
    if isempty(elements_cache{p_ID})
        [elements_cache{p_ID}, rates_cache{p_ID}] = AA279j2000_planetary_elements(p_ID);
    end

    % Turn MJD into J2000 epoch
    t_J2000 = MJD - 51544.5; % days, J2000 epoch
    t_J2000_c = t_J2000/36525; % centuries, J2000 epoch

    mu_sun = 1.3271244004193938e11; % km^3/s^2, sun

    % Determine current orbital elements based on MJD
    elements = elements_cache{p_ID} + rates_cache{p_ID}*t_J2000_c;

    % Determine position and velocity based on current orbital elements
    a = elements(1) * 149597870.7;         %   semimajor axis [AU -> km]
    e = elements(2);        %    eccentricity
    i = elements(3) ;        %   inclination [deg]
    OMEGA = elements(4) ;   %    longitude of the ascending node [deg]
    omega_bar = elements(5) ;  %  longitude of perihelion [deg]
    L = elements(6) ;  % mean longitude [deg]

    M = wrapTo360(L - omega_bar); % mean anomaly [deg]
    w = omega_bar - OMEGA; % argument of periapsis [deg]

    % Convert orbital elements to pos/vel
    [x_hc, v_hc] = Kep_to_ECI(a, e, deg2rad(i), deg2rad(w), deg2rad(OMEGA),...
                        deg2rad(M), 0, 0, mu_sun); % need radians
end

function [x, v]  = Kep_to_ECI(a, e, i, w, W, M0, t0, t, mu)
% From 279A scripts, units in radians
% inputs:
% a - semimajor axis
% e - eccentricity
% i - inclination
% w - argument of pericenter
% W - RAAN
% M0 - mean anomaly at time t0
% t0 - initial time
% t  - time
% mu  - gravitational constant (G*M)
% output: x and v in ECI inertial coordinates

% note that these elements can be singular

% number of times
n_times = numel(t);

% mean motion
n = sqrt(mu / (a .^ 3));

% mean anomaly
M = M0 + n * (t - t0);

% find eccentric anomaly
% first, let's define initial guess
if e > 0.9
    E   = pi;
else
    E = M + e * sin(M);
end

tol = 1.0e-14; % tolerance
d   = tol + 1;
iter = 0;
max_iter = 50; % guard against non-convergence at machine precision

% Newton's iterations
while d > tol && iter < max_iter
    E_old = E;
    E     = E - (E - e .* sin(E) - M) ./ (1 - e .* cos(E));
    d     = abs(E - E_old);
    iter  = iter + 1;
end

% orbit equation using eccentric anomaly
sinE = sin(E);
cosE = cos(E);

r = a .* (1 - e .* cosE);

% perifocal coordinates and velocities
r_pqw = zeros(3, n_times);
v_pqw = zeros(3, n_times);

r_pqw(1, :) = a .* (cosE - e);
r_pqw(2, :) = a .* sqrt(1 - e.^2) .* sinE;

% prefactor in velocity
v_factor = sqrt(mu .* a) ./ r;

% velocity
v_pqw(1, :) = -v_factor .* sinE;
v_pqw(2, :) =  v_factor .* sqrt(1 - e.^2) .* cosE;

% rotate to the inertial frame
R = rot(-W, 3) * rot(-i, 1) * rot(-w, 3); % compound rotation matrix

x = R * r_pqw; % rotate positions
v = R * v_pqw; % rotate velocities

end

function F = fun(r_p, mu, v_in,v_out,del_tot)
    F = asin(1/(1 + r_p*norm(v_in)^2/mu)) + ...
        asin(1/(1 + r_p*norm(v_out)^2/mu)) - del_tot;
end