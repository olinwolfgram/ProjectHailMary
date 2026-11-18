function [tIn,stateIn,tOut,stateOut,plane,burnVector,errors] = ...
    venus_flyby_arcs(vIn,vOut,rp,signedBurn,mu)

%VENUS_FLYBY_ARCS Build both hyperbolas in Venus-centered ecliptic axes.
% States are [x y z vx vy vz] in km and km/s. Time is from periapsis.

[~,eIn,aIn,turnIn,psiIn] = hyperbola_at_periapsis(norm(vIn),mu,rp);
[~,eOut,aOut,~,psiOut] = hyperbola_at_periapsis(norm(vOut),mu,rp);

n = cross(vIn,vOut);
n = n/norm(n);
u = vIn/norm(vIn);
halfTurn = turnIn/2;

% Rotate the incoming asymptote to get the tangent at periapsis.
q = cosd(halfTurn)*u+sind(halfTurn)*cross(n,u);
p = cross(q,n);
plane = [p q];
burnVector = signedBurn*q;

boundary = 3*rp; % Display boundary, not a finite SOI patch
[tIn,stateIn] = branch(eIn,aIn,rp,mu,plane,boundary,true);
[tOut,stateOut] = branch(eOut,aOut,rp,mu,plane,boundary,false);

% Check that the analytic velocities tend to the required asymptotes.
vin = velocity(-psiIn,eIn,rp,mu,plane);
vout = velocity(psiOut,eOut,rp,mu,plane);
errors = [norm(vin.'-vIn),norm(vout.'-vOut)];
assert(all(errors<1e-8),'The local hyperbolas do not match the required asymptotes.');
end

function [t,state] = branch(e,a,rp,mu,plane,boundary,isIncoming)
latus = rp*(1+e);
limit = acosd((latus/boundary-1)/e);
if isIncoming
    nu = linspace(-limit,0,201).';
else
    nu = linspace(0,limit,201).';
end

r = latus./(1+e*cosd(nu));
xyz = (r.*[cosd(nu),sind(nu)])*plane.';

% Hyperbolic Kepler equation gives the clock relative to periapsis.
H = 2*atanh(sqrt((e-1)/(e+1))*tand(nu/2));
t = sqrt((-a)^3/mu)*(e*sinh(H)-H);
state = [xyz,velocity(nu,e,rp,mu,plane)];
end

function v = velocity(nu,e,rp,mu,plane)
latus = rp*(1+e);
v = sqrt(mu/latus)*[-sind(nu),e+cosd(nu)]*plane.';
end
