function [jd,state,missKm,velocityError] = propagate_solar_leg(depJD,arrJD,r1,v1,r2,v2)
%PROPAGATE_SOLAR_LEG Integrate Sun-only motion and check the Lambert endpoint.

% Dates are Julian days; positions are in km and velocities are in km/s.
mu = body_constants('Sun');
dt = (arrJD-depJD)*86400;

assert(dt>0,'Solar leg arrival must follow departure.');

times = linspace(0,dt,ceil(dt/8640)+1); % Output every 0.1 day or less

opts = odeset('RelTol',1e-11,'AbsTol',1e-9);
[t,state] = ode113(@(~,x)[x(4:6);-mu*x(1:3)/norm(x(1:3))^3],times,[r1;v1],opts);

jd = depJD+t/86400;

missKm = norm(state(end,1:3).'-r2);
velocityError = norm(state(end,4:6).'-v2);
assert(missKm<10 && velocityError<1e-5, ...
    'Solar leg did not close: %.3g km position error, %.3g km/s velocity error.', ...
    missKm,velocityError);
end
