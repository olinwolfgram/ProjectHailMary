function [rp,signedBurn,turnDeg] = powered_swingby(vIn,vOut,mu)
%POWERED_SWINGBY Solve for the shared vacuum periapsis and its tangential burn.
% Lecture 9, pp. 29-30. Inputs are km/s and km^3/s^2; rp is in km.

vIn = vIn(:);
vOut = vOut(:);
assert(numel(vIn)==3 && numel(vOut)==3 && isreal([vIn;vOut]) && all(isfinite([vIn;vOut])), ...
    'Use two finite, real three-component excess velocities.');
assert(isfinite(mu) && mu>0,'Use a positive gravity parameter.');

speedIn = norm(vIn);
speedOut = norm(vOut);
assert(speedIn>0 && speedOut>0,'Both excess speeds must be nonzero.');
turn = atan2(norm(cross(vIn,vOut)),dot(vIn,vOut));
assert(turn>0 && turn<pi,'A finite periapsis requires a turn between 0 and 180 deg.');

% Add the two half-turns and solve for the radius that matches the vectors.
halfTurn = @(r,s) asin(1./(1+r*s^2/mu));
residual = @(r) halfTurn(r,speedIn)+halfTurn(r,speedOut)-turn;
upper = 2*(1/sin(turn/2)-1)*mu/min(speedIn,speedOut)^2;
rp = fzero(residual,[0 upper],optimset('TolX',1e-12));

vpIn = hyperbola_at_periapsis(speedIn,mu,rp);
vpOut = hyperbola_at_periapsis(speedOut,mu,rp);
signedBurn = vpOut-vpIn;
turnDeg = rad2deg(turn);
end
