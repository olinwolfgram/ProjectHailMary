function mission = run_venus_gram_flyby(options)
%RUN_VENUS_GRAM_FLYBY Compare atmospheric and vacuum Earth-Venus-Earth legs.
% Run without inputs for the 111 km example from Pearl's low-flyby search.
% Supply three UTC Julian dates or one row from Pearl's candidate table.

arguments
    options.DatesJD (1,:) double {mustBeFinite} = []
    options.PearlCandidate table = table()
    options.GramRoot (1,1) string = ""
    options.Beta (1,1) double {mustBePositive,mustBeFinite} = 100 % kg/m^2
    options.ParkingAltitudeKm (1,1) double {mustBeNonnegative,mustBeFinite} = 200
    options.InterfaceKm (1,1) double {mustBePositive,mustBeFinite} = 300
    options.MinAltitudeKm (1,1) double {mustBePositive,mustBeFinite} = 80
    options.ProfileStepKm (1,1) double {mustBePositive,mustBeFinite} = 0.5
    options.PatchRadiusKm (1,1) double {mustBePositive,mustBeFinite} = 616000
    options.MaxStepSeconds (1,1) double {mustBePositive,mustBeFinite} = 5
    options.Verbose (1,1) logical = true
end

root = fileparts(mfilename('fullpath'));
pearl = fullfile(fileparts(root),"Pearl's Code");
assert(isfolder(pearl),"Keep venus_gram alongside Pearl's Code in the repository.");
addpath(pearl);

gramRoot = options.GramRoot;
if gramRoot==""
    gramRoot = string(getenv('VENUS_GRAM_ROOT'));
end
if gramRoot==""
    gramRoot = fullfile(root,'runtime');
end

dates = options.DatesJD;
candidate = options.PearlCandidate;
if ~isempty(candidate)
    assert(isempty(dates),'Supply either PearlCandidate or DatesJD.');
    columns = {'dep_MJD','fb_MJD','ret_MJD'};
    assert(height(candidate)==1 && all(ismember(columns,candidate.Properties.VariableNames)), ...
        'PearlCandidate must be one table row with dep_MJD, fb_MJD and ret_MJD.');
    dates = candidate{1,columns}+2400000.5;
elseif isempty(dates)
    dates = juliandate(datetime([2024 12 24;2025 4 25;2026 1 6],'TimeZone','UTC')).';
end
assert(isnumeric(dates) && isrow(dates) && numel(dates)==3 && ...
    all(isfinite(dates)) && all(diff(dates)>0), ...
    'DatesJD must contain Earth departure, Venus periapsis and Earth arrival in order.');
assert(options.MinAltitudeKm < options.InterfaceKm, ...
    'The altitude floor must be below the atmosphere interface.');

[mu, radius] = body_constants('Venus');
muSun = body_constants('Sun');
[rVenus, vVenus] = planet_state(2,dates(2));
[rLaunch, vEarthLaunch] = planet_state(3,dates(1));
[rEarth, vEarth] = planet_state(3,dates(3));
[vLaunch, vArrival, solved] = pearlLambert(muSun,rLaunch,rVenus,86400*(dates(2)-dates(1)));
assert(solved,'Earth-to-Venus Lambert solve failed.');
[vDepart, vEarthArrival, solved] = pearlLambert(muSun,rVenus,rEarth,86400*(dates(3)-dates(2)));
assert(solved,'Vacuum Venus-to-Earth Lambert solve failed.');
injection = earth_injection_delta_v(norm(vLaunch-vEarthLaunch)^2,options.ParkingAltitudeKm);
vIn = vArrival-vVenus;
target = vDepart-vVenus;
[rp, signedBurn] = powered_swingby(vIn,target,mu);
[~,~,~,~,plane] = venus_flyby_arcs(vIn,target,rp,signedBurn,mu);
assert(rp>radius+options.MinAltitudeKm && rp<radius+options.InterfaceKm, ...
    'Vacuum periapsis must lie between the altitude floor and the interface.');
assert(options.PatchRadiusKm>radius+options.InterfaceKm, ...
    'The return patch radius must be above the interface.');

[~,e,a] = hyperbola_at_periapsis(norm(vIn),mu,rp);
p = rp*(1+e);
entryRadius = radius+options.InterfaceKm;
nu = -acos((p/entryRadius-1)/e);
H = 2*atanh(sqrt((e-1)/(e+1))*tan(nu/2));
tEntry = sqrt((-a)^3/mu)*(e*sinh(H)-H);
xEntry = [entryRadius*plane*[cos(nu);sin(nu)]; ...
    sqrt(mu/p)*plane*[-sin(nu);e+cos(nu)]];

% IAU_VENUS rotation from the bundled pck00011.tpc. The UTC-to-TDB
% approximation here only sets the profile longitude, not the ODE clock.
[toFixed, pole, spin] = venusOrientation(dates(2));
fixed = toFixed*plane(:,1);
latitude = atan2d(fixed(3),hypot(fixed(1),fixed(2)));
longitude = mod(atan2d(fixed(2),fixed(1)),360);
profile = readGramProfile(gramRoot,options,dates(2),latitude,longitude);
rhoAt = griddedInterpolant(profile.Height_km,log(profile.Density_kgm3),'linear','nearest');
eastAt = griddedInterpolant(profile.Height_km,profile.EWWind_ms,'linear','nearest');
northAt = griddedInterpolant(profile.Height_km,profile.NSWind_ms,'linear','nearest');

base = odeset('RelTol',1e-10,'AbsTol',1e-11,'MaxStep',options.MaxStepSeconds);
mission.DatesJD = dates;
mission.Beta_kgm2 = options.Beta;
mission.ParkingAltitude_km = options.ParkingAltitudeKm;
mission.EarthInjection_kms = injection;
mission.VenusRadius_km = radius;
mission.AtmosphereBand_km = [options.MinAltitudeKm options.InterfaceKm];
mission.FlybyPlane = plane;
mission.VacuumPeriapsis_km = rp-radius;
mission.TargetVinf_kms = target;
mission.ProfileLatitude_deg = latitude;
mission.ProfileLongitudeEast_deg = longitude;
mission.Profile = profile;
[mission.Outbound.JD,mission.Outbound.State] = propagate_solar_leg( ...
    dates(1),dates(2),rLaunch,vLaunch,rVenus,vArrival);
[mission.VacuumReturn.JD,mission.VacuumReturn.State] = propagate_solar_leg( ...
    dates(2),dates(3),rVenus,vDepart,rEarth,vEarthArrival);
mission.VacuumBallistic = propagatePass(false);
mission.DragBallistic = propagatePass(true);
mission.VacuumPatch = patchReturn(mission.VacuumBallistic);
mission.DragPatch = patchReturn(mission.DragBallistic);
if options.Verbose
    printResults(mission);
end

    function [derivative,q,dragNorm] = rhs(~,x,hasDrag)
        r = x(1:3);
        v = x(4:6);
        distance = norm(r);
        drag = zeros(3,1);
        q = 0;
        if hasDrag
            altitude = distance-radius;
            up = r/distance;
            east = cross(pole,up);
            if norm(east)<1e-12
                east = toFixed(2,:).';
            else
                east = east/norm(east);
            end
            north = cross(up,east);
            winds = (eastAt(altitude)*east+northAt(altitude)*north)/1000;
            rel = v-spin*cross(pole,r)-winds;
            q = 0.5*exp(rhoAt(altitude))*(1000*norm(rel))^2;
            drag = -q/options.Beta/1000 * rel/max(norm(rel),realmin);
        end
        derivative = [v; -mu*r/distance^3+drag];
        dragNorm = norm(drag);
    end

    function [value,terminal,direction] = entryEvents(~,x)
        value = [dot(x(1:3),x(4:6)); norm(x(1:3))-radius-options.MinAltitudeKm];
        terminal = [1;1];
        direction = [1;-1];
    end

    function [value,terminal,direction] = exitEvents(~,x)
        distance = norm(x(1:3));
        value = [distance-entryRadius; distance-radius-options.MinAltitudeKm];
        terminal = [1;1];
        direction = [1;-1];
    end

    function arc = propagatePass(hasDrag)
        dynamics = @(t,x)rhs(t,x,hasDrag);
        [arc.Time_s, arc.State, ~, ~, inboundEvent] = ode113(dynamics, ...
            [tEntry 14400],xEntry,odeset(base,'Events',@entryEvents));
        arc.Status = "no periapsis";
        arc.ExitVinf_kms = nan(3,1);
        if any(inboundEvent==2)
            arc.Status = "altitude floor";
        elseif any(inboundEvent==1)
            [t,x,~,~,event] = ode113(dynamics,[arc.Time_s(end) 14400], ...
                arc.State(end,:).',odeset(base,'Events',@exitEvents));
            arc.Time_s = [arc.Time_s;t(2:end)];
            arc.State = [arc.State;x(2:end,:)];
            arc.Status = "no outbound crossing";
            if any(event==2)
                arc.Status = "altitude floor";
            elseif any(event==1)
                arc.ExitVinf_kms = outgoingVinf(x(end,:).',mu);
                arc.Status = "bound at exit";
                if all(isfinite(arc.ExitVinf_kms))
                    arc.Status = "escaped";
                end
            end
        end
        arc.MinAltitude_km = min(vecnorm(arc.State(:,1:3),2,2))-radius;
        speed = vecnorm(arc.State(:,4:6),2,2);
        arc.ExitEnergy_km2s2 = speed(end)^2/2-mu/norm(arc.State(end,1:3));
        arc.DynamicPressure_Pa = zeros(size(arc.Time_s));
        arc.Drag_ms2 = zeros(size(arc.Time_s));
        if hasDrag
            for k = 1:size(arc.State,1)
                [~,q,dragNorm] = rhs(arc.Time_s(k),arc.State(k,:).',true);
                arc.DynamicPressure_Pa(k) = q;
                arc.Drag_ms2(k) = 1000*dragNorm;
            end
        end
        arc.PeakDynamicPressure_Pa = max(arc.DynamicPressure_Pa);
        arc.PeakDrag_ms2 = max(arc.Drag_ms2);
    end

    function patch = patchReturn(arc)
        patch = struct('Status',arc.Status,'Burn_kms',nan(3,1), ...
            'Radius_km',NaN,'JD',NaN,'EarthArrivalVinf_kms',NaN,'EarthMiss_km',NaN, ...
            'CoastTime_s',zeros(0,1),'CoastState',zeros(0,6), ...
            'ReturnJD',zeros(0,1),'ReturnState',zeros(0,6));
        if arc.Status~="escaped"
            return
        end
        [patch.CoastTime_s,patch.CoastState,te,xe] = ode113(@(t,x)rhs(t,x,false), ...
            [arc.Time_s(end) arc.Time_s(end)+7*86400],arc.State(end,:).', ...
            odeset(base,'MaxStep',1200,'Events',@patchEvent));
        patch.Status = "no finite patch crossing";
        if isempty(te)
            return
        end
        patch.JD = dates(2)+te(end)/86400;
        patch.Radius_km = norm(xe(1:3));
        [rv,vv] = planet_state(2,patch.JD);
        r0 = rv+xe(1:3).';
        v0 = vv+xe(4:6).';
        [v1,v2,solvedReturn] = pearlLambert(muSun,r0,rEarth,(dates(3)-patch.JD)*86400);
        patch.Status = "return Lambert failed";
        if ~solvedReturn
            return
        end
        patch.Status = "targeted";
        patch.Burn_kms = v1-v0;
        patch.EarthArrivalVinf_kms = norm(v2-vEarth);
        [patch.ReturnJD,patch.ReturnState,patch.EarthMiss_km] = ...
            propagate_solar_leg(patch.JD,dates(3),r0,v1,rEarth,v2);
    end

    function [value,terminal,direction] = patchEvent(~,x)
        value = norm(x(1:3))-options.PatchRadiusKm;
        terminal = 1;
        direction = 1;
    end

end

function [v1,v2,solved] = pearlLambert(mu,r1,r2,seconds)
v1 = nan(3,1);
v2 = nan(3,1);
solved = false;
if seconds>0 && norm(cross(r1,r2))>1e-8*norm(r1)*norm(r2)
    [v1,v2,status] = AA279lambert_curtis(mu,r1,r2,'pro',0,seconds);
    solved = strcmp(status,'finished') && isreal([v1;v2]) && all(isfinite([v1;v2]));
end
end

function vInf = outgoingVinf(x,mu)
r = x(1:3);
v = x(4:6);
energy = dot(v,v)/2-mu/norm(r);
vInf = nan(3,1);
if energy>0
    h = cross(r,v);
    eVector = cross(v,h)/mu-r/norm(r);
    e = norm(eVector);
    p = eVector/e;
    q = cross(h/norm(h),p);
    vInf = sqrt(2*energy)*(-p/e+sqrt(1-1/e^2)*q);
end
end

function [toFixed,pole,spin] = venusOrientation(jd)
% The LSK in this GRAM distribution gives TT-UTC=69.184 s for the
% project epochs. TDB-TT is under 2 ms and omitted for this profile query.
days = jd-2451545+69.184/86400;
spin = deg2rad(-1.4813688)/86400;
rotZ = @(angle)[cosd(angle) sind(angle) 0; -sind(angle) cosd(angle) 0; 0 0 1];
rotX = @(angle)[1 0 0; 0 cosd(angle) sind(angle); 0 -sind(angle) cosd(angle)];
toFixed = rotZ(160.20-1.4813688*days)*rotX(90-67.16)* ...
    rotZ(90+272.76)*rotX(-23.439291111);
pole = toFixed(3,:).';
end

function profile = readGramProfile(gramRoot,options,jd,latitude,longitude)
[found, location] = fileattrib(gramRoot);
assert(found && location.directory, ...
    'GRAM runtime missing at %s. Set GramRoot or the VENUS_GRAM_ROOT environment variable.',gramRoot);
exe = fullfile(location.Name,'Windows','VenusGRAM.exe');
spice = fullfile(location.Name,'SPICE');
assert(isfile(exe) && isfolder(spice), ...
    'Set GramRoot to the GRAM directory containing Windows and SPICE.');
assert(~contains(spice,["'" newline]),'The GRAM data path contains a namelist delimiter.');

scratch = tempname;
mkdir(scratch);
original = pwd;
cleanup = onCleanup(@()removeScratch(original,scratch)); %#ok<NASGU>
cd(scratch);
start = datetime(jd,'ConvertFrom','juliandate','TimeZone','UTC');
[y,m,d] = ymd(start);
[hr,mn,sec] = hms(start);
count = ceil((options.InterfaceKm-options.MinAltitudeKm)/options.ProfileStepKm)+1;
input = fopen('profile.txt','wt');
assert(input>=0,'Could not create the temporary Venus-GRAM namelist.');
closeInput = onCleanup(@()fclose(input));
fprintf(input,[' $INPUT\n SpicePath = ''%s''\n ListFileName = ''profile_LIST''\n', ...
    ' ColumnFileName = ''profile_OUTPUT''\n EastLongitudePositive = 1\n', ...
    ' IsPlanetoCentric = 1\n TimeFrame = 0\n TimeScale = 1\n', ...
    ' Year = %d\n Month = %d\n Day = %d\n Hour = %d\n Minute = %d\n Seconds = %.7f\n', ...
    ' NumberOfPositions = %d\n InitialHeight = %.9f\n DeltaHeight = %.9f\n', ...
    ' InitialLatitude = %.9f\n InitialLongitude = %.9f\n', ...
    ' DensityPerturbationScale = 0\n EWWindPerturbationScale = 0\n', ...
    ' NSWindPerturbationScale = 0\n NumberOfMonteCarloRuns = 1\n', ...
    ' ExtraPrecision = 1\n UseLegacyOutputs = 0\n $END\n'], ...
    spice,y,m,d,hr,mn,sec,count,options.MinAltitudeKm,options.ProfileStepKm,latitude,longitude);
clear closeInput
[status,output] = system(sprintf('"%s" -file profile.txt',exe));
assert(status==0 && isfile('profile_OUTPUT.csv'),'Venus-GRAM failed:\n%s',output);
profile = readtable('profile_OUTPUT.csv');
assert(height(profile)==count && all(diff(profile.Height_km)>0) && ...
    all(isfinite(profile.Density_kgm3) & profile.Density_kgm3>0) && ...
    all(isfinite(profile.EWWind_ms) & isfinite(profile.NSWind_ms)), ...
    'Venus-GRAM returned an incomplete or invalid atmosphere profile.');
end

function removeScratch(original,scratch)
cd(original);
rmdir(scratch,'s');
end

function printResults(mission)
fprintf('Venus-GRAM flyby at beta %.1f kg/m^2\n',mission.Beta_kgm2);
fprintf('Earth injection from %.0f km parking orbit: %.3f km/s\n', ...
    mission.ParkingAltitude_km,mission.EarthInjection_kms);
fprintf('Nominal vacuum periapsis %.2f km\n',mission.VacuumPeriapsis_km);
drag = mission.DragBallistic;
fprintf('Drag-only pass: %s | minimum altitude %.2f km | peak q %.3g Pa | peak drag %.3g m/s^2\n', ...
    drag.Status,drag.MinAltitude_km,drag.PeakDynamicPressure_Pa,drag.PeakDrag_ms2);
fprintf('Vacuum comparison at the finite patch: %s | correction %.2f m/s\n', ...
    mission.VacuumPatch.Status,1000*norm(mission.VacuumPatch.Burn_kms));
fprintf('Post-flyby correction with drag: %s | burn %.2f m/s | Earth arrival v_inf %.3f km/s\n', ...
    mission.DragPatch.Status,1000*norm(mission.DragPatch.Burn_kms), ...
    mission.DragPatch.EarthArrivalVinf_kms);
if mission.DragPatch.Status=="targeted"
    fprintf('Post-flyby burn %.2f hours after nominal periapsis | injection + correction %.3f km/s\n', ...
        mission.DragPatch.CoastTime_s(end)/3600,mission.EarthInjection_kms+norm(mission.DragPatch.Burn_kms));
end
end
