function [deltaV,rp] = earth_injection_delta_v(c3,parkingAltitudeKm)
%EARTH_INJECTION_DELTA_V Tangential escape burn from a circular parking orbit.

[mu,radius] = body_constants('Earth');
rp = radius+parkingAltitudeKm;

deltaV = sqrt(c3+2*mu/rp)-sqrt(mu/rp);
end
