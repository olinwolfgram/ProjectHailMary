clear; close all; clc;

% Constants

% Grav parameter (km^3/s^2)
mu.venus = 324858.630;
mu.earth = 398600.440;
mu.sun = 1.3271244004193938e11;

% Mean radius (km)
r_m.venus = 6051.80;
r_m.earth = 6371.01;

% Arrival and Departures
p1_ID = 3; % Earth
p2_ID = 2; % Venus
dep_date = [2020, 01, 01, 00, 00, 00]; % Y, M, D, H, M, S

tStart = mjuliandate(datetime(dep_date, 'TimeZone','UTC'));
step = 10; % days, so my computer doesn't explode rn

% Modified Julian Dates used
dep_dates = tStart : step : tStart + 365*5;
arr_times = 15: step : 10*30; % 15 days to 8 months

% Iterate over MJDs chosen
p1_C3_total_grid = zeros(length(dep_dates),length(arr_times));
p2_C3_total_grid = p1_C3_total_grid;
DLA1_total_grid = p1_C3_total_grid;
DLA2_total_grid = p1_C3_total_grid;
for j = 1:length(dep_dates) % Iterate over departures
    dep_MJD = dep_dates(j);
    for k = 1:length(arr_times) % Iterate over arrivals
        arr_MJD = dep_MJD + arr_times(k);
        [p1_Vinf, p2_Vinf] = p1_To_p2_dV(p1_ID, p2_ID, dep_MJD, arr_MJD);

        % Calculate Vinf (CHANGE)
        p1_C3_total_grid(j,k) = norm(p1_Vinf)^2;
        p2_C3_total_grid(j,k) = norm(p2_Vinf)^2;

        % Calculate launch asymptote
        DLA1_total_grid(j,k) = asind(p1_Vinf(3)/norm(p1_Vinf));
        DLA2_total_grid(j,k) = asind(p2_Vinf(3)/norm(p2_Vinf));

    end
end

% Grid of departure MJD and time of flight (rows = TOF, cols = departure)
[DEP, TOF] = meshgrid(dep_dates, arr_times);
ARR = DEP + TOF;                 % arrival MJD at every grid point

% Convert MJD -> MATLAB datenum so the axes can show calendar dates
mjd2datenum = 678942;            % datenum of MJD 0 (17-Nov-1858)
DEP_dn = DEP + mjd2datenum;
ARR_dn = ARR + mjd2datenum;

% Contour levels: focus on the useful low-dV region
p1_C3_min = min(p1_C3_total_grid(:))
p1_levels = ceil(p1_C3_min) : 2 : p1_C3_min + 12; % km^2/s^2

figure;
[C, h] = contour(DEP_dn, ARR_dn, p1_C3_total_grid.', p1_levels, 'ShowText', 'on');
colorbar; clim([p1_C3_min, p1_C3_min + 10]);
datetick('x', 'mmm-yy', 'keeplimits');
datetick('y', 'mmm-yy', 'keeplimits');
xlabel('Departure Date');
ylabel('Arrival Date');
title('Earth \rightarrow Venus Departure C3 (km^2/s^2)');
grid on;

p2_C3_min = min(p2_C3_total_grid(:))
p2_levels = ceil(p2_C3_min) : 2 : p2_C3_min + 12;   % km^2/s^2

figure;
[C, h] = contour(DEP_dn, ARR_dn, p2_C3_total_grid.', p2_levels, 'ShowText', 'on');
colorbar; clim([p2_C3_min, p2_C3_min + 10]);
datetick('x', 'mmm-yy', 'keeplimits');
datetick('y', 'mmm-yy', 'keeplimits');
xlabel('Departure Date');
ylabel('Arrival Date');
title('Earth \rightarrow Venus Arrival C3 (km^2/s^2)');
grid on;

p12_C3_min = min(p1_C3_total_grid(:).' + p2_C3_total_grid(:).')
p12_levels = ceil(p12_C3_min) : 2 : ceil(p12_C3_min) + 12;   % km^2/s^2

figure;
[C, h] = contour(DEP_dn, ARR_dn, p1_C3_total_grid.' + p2_C3_total_grid.', p12_levels, 'ShowText', 'on');
colorbar; clim([p12_C3_min, p12_C3_min + 12]);
datetick('x', 'mmm-yy', 'keeplimits');
datetick('y', 'mmm-yy', 'keeplimits');
xlabel('Departure Date');
ylabel('Arrival Date');
title('Earth \rightarrow Venus C3 Total (km^2/s^2)');
grid on;


%% Functions -------------------
% Interplanetary Lambert Wrapper
function [p1_Vinf, p2_Vnif] = p1_To_p2_dV(p1_ID, p2_ID, dep_MJD, arr_MJD)

% Determine location of Earth at departure and Venus at arrival
[p1.x_dep_hc, p1.v_dep_hc] = heliocentric(p1_ID, dep_MJD);
[p2.x_arr_hc, p2.v_arr_hc] = heliocentric(p2_ID, arr_MJD);

% Solve Lambert
dm = 'pro';
nrev = 0;
dt = arr_MJD - dep_MJD; % days
mu_sun = 1.3271244004193938e11; % km^3/s^2, sun
[v1_out, v2_out, error_out] = AA279lambert_curtis(mu_sun,...
                            p1.x_dep_hc, p2.x_arr_hc, dm, nrev, dt*24*3600);

% Find planetocentric Vinf
p1_Vinf = v1_out - p1.v_dep_hc;
p2_Vnif = v2_out - p2.v_arr_hc;
end

function [x_hc, v_hc] = heliocentric(p_ID, MJD)
    % Turn MJD into J2000 epoch
    t_J2000 = MJD - 51544.5; % days, J2000 epoch
    t_J2000_c = t_J2000/36525; % centuries, J2000 epoch

    mu_sun = 1.3271244004193938e11; % km^3/s^2, sun

    [elements, rates] = AA279j2000_planetary_elements(p_ID); 

    % Determine current orbital elements based on MJD
    elements = elements + rates*t_J2000_c;

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

% Newton's iterations
while d > tol
    E_old = E;
    E     = E - (E - e .* sin(E) - M) ./ (1 - e .* cos(E));
    d     = abs(E - E_old);
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