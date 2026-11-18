function [mu,radius,au] = body_constants(body)
%BODY_CONSTANTS Body values in km and seconds.

au = 149597870.7; % km per AU

switch body
    case 'Sun'
        mu = 1.3271244004193938e11;
        radius = []; % not fr using this, just empty for now
    case 'Earth'
        mu = 398600.440;
        radius = 6378.0;
    case 'Venus'
        mu = 324858.630;
        radius = 6051.80;
    otherwise
        error('Choose Sun, Earth, or Venus.');
end
end
