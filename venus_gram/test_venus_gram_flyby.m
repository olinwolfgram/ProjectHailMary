function tests = test_venus_gram_flyby
%TEST_VENUS_GRAM_FLYBY End-to-end checks using the bundled Windows GRAM exe.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
candidate = table(60668,60790,61046, ...
    'VariableNames',{'dep_MJD','fb_MJD','ret_MJD'});
testCase.TestData.Mission = run_venus_gram_flyby(PearlCandidate=candidate,Verbose=false);
end

function testPearlHandoff(testCase)
m = testCase.TestData.Mission;
verifyEqual(testCase,m.DatesJD,[2460668.5 2460790.5 2461046.5]);
verifyEqual(testCase,m.VacuumPeriapsis_km,111.11,'AbsTol',0.02);
parent = fileparts(fileparts(mfilename('fullpath')));
pearl = fullfile(parent,"Pearl's Code");
for helper = ["AA279lambert_curtis","AA279j2000_planetary_elements","rot"]
    verifyEqual(testCase,string(which(helper)),fullfile(pearl,helper+".m"));
end
end

function testVacuumClosure(testCase)
m = testCase.TestData.Mission;
vacuum = m.VacuumBallistic;
mu = body_constants('Venus');
verifyEqual(testCase,vacuum.Status,"escaped");
verifyEqual(testCase,vacuum.MinAltitude_km,m.VacuumPeriapsis_km,'AbsTol',1e-4);
energy = sum(vacuum.State(:,4:6).^2,2)/2-mu./vecnorm(vacuum.State(:,1:3),2,2);
verifyEqual(testCase,energy,energy(1)*ones(size(energy)),'AbsTol',1e-7);
verifyEqual(testCase,m.VacuumBallistic.PeakDynamicPressure_Pa,0);
verifyEqual(testCase,m.VacuumBallistic.PeakDrag_ms2,0);
end

function testAtmosphericPass(testCase)
m = testCase.TestData.Mission;
verifyGreaterThan(testCase,height(m.Profile),100);
verifyTrue(testCase,all(m.Profile.Density_kgm3>0));
verifyEqual(testCase,m.DragBallistic.Status,"escaped");
verifyGreaterThan(testCase,m.DragBallistic.ExitEnergy_km2s2,0);
verifyLessThan(testCase,m.DragBallistic.ExitEnergy_km2s2,m.VacuumBallistic.ExitEnergy_km2s2);
verifyGreaterThan(testCase,m.DragBallistic.PeakDynamicPressure_Pa,0);
verifyEqual(testCase,m.DragBallistic.PeakDrag_ms2, ...
    m.DragBallistic.PeakDynamicPressure_Pa/m.Beta_kgm2,'RelTol',1e-12);
verifySize(testCase,m.DragBallistic.DynamicPressure_Pa,size(m.DragBallistic.Time_s));
verifyEqual(testCase,max(m.DragBallistic.DynamicPressure_Pa),m.DragBallistic.PeakDynamicPressure_Pa);
end

function testUnpoweredCoast(testCase)
m = testCase.TestData.Mission;
verifyEqual(testCase,m.DragBallistic.State(1,:),m.VacuumBallistic.State(1,:));
verifyTrue(testCase,all(diff(m.DragBallistic.Time_s)>0));
verifyEqual(testCase,m.DragBallistic.State(end,:),m.DragPatch.CoastState(1,:));
verifyEqual(testCase,m.DragBallistic.Time_s(end),m.DragPatch.CoastTime_s(1));
verifyGreaterThan(testCase,m.DragPatch.CoastTime_s(end),m.DragBallistic.Time_s(end));
end

function testFinitePatchReturn(testCase)
m = testCase.TestData.Mission;
verifyEqual(testCase,m.VacuumPatch.Status,"targeted");
verifyEqual(testCase,m.DragPatch.Status,"targeted");
verifyEqual(testCase,m.DragPatch.Radius_km,616000,'AbsTol',1e-3);
verifyGreaterThan(testCase,m.DragPatch.JD,m.DatesJD(2));
verifyLessThan(testCase,m.DragPatch.JD,m.DatesJD(3));
verifyLessThan(testCase,m.DragPatch.EarthMiss_km,1);
verifyGreaterThan(testCase,norm(m.DragPatch.Burn_kms-m.VacuumPatch.Burn_kms),1e-4);
[rVenus,vVenus] = planet_state(2,m.DragPatch.JD);
coastEnd = m.DragPatch.CoastState(end,:).';
returnStart = m.DragPatch.ReturnState(1,:).';
verifyEqual(testCase,returnStart(1:3),rVenus+coastEnd(1:3),'AbsTol',1e-6);
verifyEqual(testCase,returnStart(4:6),vVenus+coastEnd(4:6)+m.DragPatch.Burn_kms,'AbsTol',1e-10);
end

function testConvergence(testCase)
m = testCase.TestData.Mission;
refined = run_venus_gram_flyby(DatesJD=m.DatesJD,InterfaceKm=350, ...
    ProfileStepKm=0.25,MaxStepSeconds=2.5,Verbose=false);
verifyEqual(testCase,refined.DragPatch.Status,"targeted");
verifyLessThan(testCase,norm(m.DragPatch.Burn_kms-refined.DragPatch.Burn_kms),1e-3);
verifyEqual(testCase,m.DragBallistic.MinAltitude_km,refined.DragBallistic.MinAltitude_km,'AbsTol',0.01);
end

function testAltitudeFloor(testCase)
m = run_venus_gram_flyby(Beta=1,MinAltitudeKm=110.9,Verbose=false);
verifyEqual(testCase,m.DragBallistic.Status,"altitude floor");
verifyEqual(testCase,m.DragBallistic.MinAltitude_km,110.9,'AbsTol',1e-4);
verifyEqual(testCase,m.DragPatch.Status,"altitude floor");
verifyTrue(testCase,all(isnan(m.DragPatch.Burn_kms)));
end

function testMissionPlots(testCase)
existing = get(groot,'DefaultFigureVisible');
testCase.addTeardown(@()set(groot,'DefaultFigureVisible',existing));
set(groot,'DefaultFigureVisible','off');
figures = plot_venus_gram_flyby(testCase.TestData.Mission);
testCase.addTeardown(@()delete(figures));
verifyNumElements(testCase,figures,2);
verifyNumElements(testCase,findobj(figures(1),'Type','axes'),5);
verifyNumElements(testCase,findobj(figures(2),'Type','axes'),2);
for view = ["VenusOverview","VenusCloseup"]
    ax = findobj(figures(1),'Type','axes','Tag',char(view));
    testCase.assertNumElements(ax,1);
    verifyEqual(testCase,ax.DataAspectRatio,[1 1 1]);
    disk = findobj(ax,'Type','patch','Tag','VenusDisk');
    testCase.assertNumElements(disk,1);
    distance = hypot(disk.XData,disk.YData);
    verifyEqual(testCase,distance, ...
        testCase.TestData.Mission.VenusRadius_km*ones(size(distance)),'AbsTol',1e-8);
end
end
