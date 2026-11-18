function [r, v] = planet_state(planet_id, JD)
% Heliocentric position [km] and velocity [km/s], as 3-by-1 vectors.
% Frame: J2000 ecliptic. planet_id: 1 = Mercury, ..., 9 = Pluto.
% JD is the Julian date, not the modified Julian date.

[mu,~,AU] = body_constants('Sun');

% Advance the J2000 elements to the requested date.
[elements, rates] = AA279j2000_planetary_elements(planet_id);
T = (JD - 2451545.0)/36525; % Julian centuries since J2000
elements = elements + rates*T;

a = elements(1)*AU;
e = elements(2);
inc = deg2rad(elements(3));
Omega = deg2rad(elements(4));
omega = deg2rad(elements(5) - elements(4));
M = deg2rad(mod(elements(6) - elements(5), 360));

% Solve M = E - e*sin(E) using Newton's method.
E = M;
for k = 1:50
    dE = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - dE;
    if abs(dE) < 1e-12
        break
    end
end

% Position and two-body velocity in the orbital plane.
r_orbit = a*[cos(E) - e; sqrt(1 - e^2)*sin(E); 0];
v_orbit = sqrt(mu/a)/(1 - e*cos(E))* ...
          [-sin(E); sqrt(1 - e^2)*cos(E); 0];

% Pearl's rot uses passive rotations, hence the negative angles.
Q = rot(-Omega,3)*rot(-inc,1)*rot(-omega,3);

r = Q*r_orbit;
v = Q*v_orbit;
end
