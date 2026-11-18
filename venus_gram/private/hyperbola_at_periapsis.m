function [vp,e,a,turn,psi,b,energy] = hyperbola_at_periapsis(vinf,mu,rp)
%HYPERBOLA_AT_PERIAPSIS Lecture 9, pp. 13, 19-24 (km and seconds).
% These inputs set the magnitudes, but not the hyperbola's orientation.

assert(all(isfinite([vinf mu rp])) && vinf>0 && mu>0 && rp>0, ...
    'Excess speed, gravity parameter, and radius must be finite and positive.');

vp = sqrt(vinf^2+2*mu/rp); % Energy conservation from infinity to periapsis
e = 1+rp*vinf^2/mu;
a = -mu/vinf^2;

turn = 2*asind(1/e);
psi = acosd(-1/e);
b = rp*vp/vinf; % b*v_infinity = r_p*v_p
energy = vinf^2/2;
end
