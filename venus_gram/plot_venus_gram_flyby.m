function figures = plot_venus_gram_flyby(mission)
%PLOT_VENUS_GRAM_FLYBY Plot the local drag pass and its E–V–E return.
vacuumColor = [0.19 0.30 0.45];
dragColor = [0.80 0.31 0.15];
bandColor = [0.85 0.93 0.97];
vacuum = mission.VacuumBallistic;
drag = mission.DragBallistic;
radius = mission.VenusRadius_km;
band = mission.AtmosphereBand_km;
plane = mission.FlybyPlane;
bandLabel = sprintf('GRAM band, %g–%g km altitude',band);
circle = [sind((0:0.5:360).'),cosd((0:0.5:360).')];
muVenus = body_constants('Venus');

% Extend the common approach with Venus gravity for the planet-scale view.
entry = vacuum.State(1,:).';
backTime = 2*radius/norm(entry(4:6));
[~,inbound] = ode113(@(~,x)[x(4:6);-muVenus*x(1:3)/norm(x(1:3))^3], ...
    linspace(0,-backTime,151),entry,odeset('RelTol',1e-10,'AbsTol',1e-11));
approach = inbound(:,1:3)*plane;
vacXY = [vacuum.State(:,1:3);mission.VacuumPatch.CoastState(:,1:3)]*plane;
dragXY = [drag.State(:,1:3);mission.DragPatch.CoastState(:,1:3)]*plane;

figures(1) = figure('Name','Venus atmospheric pass','Color','white', ...
    'Position',[65 55 1430 790]);
layout = tiledlayout(figures(1),2,3,'TileSpacing','compact','Padding','compact');
title(layout,sprintf('Venus flyby  |  %.2f km closest approach  |  Ballistic coefficient %g kg/m²', ...
    drag.MinAltitude_km,mission.Beta_kgm2));

ax = nexttile(layout,1,[2 1]);
ax.Tag = 'VenusOverview';
drawFlyby(ax);
xlim(ax,[-1.17 1.17]*radius);
ylim(ax,[-1.12 1.23]*radius);
text(ax,0,0,'Venus','HorizontalAlignment','center','FontSize',16, ...
    'FontWeight','bold','Color',[0.38 0.28 0.14]);
title(ax,'Planet-scale view');
legend(ax,'Location','southoutside','NumColumns',1);

ax = nexttile(layout,2);
ax.Tag = 'VenusCloseup';
drawFlyby(ax);
halfWidth = max(500,2*band(2));
xlim(ax,[-halfWidth halfWidth]);
ylim(ax,radius+[-0.30*halfWidth band(2)+0.10*halfWidth]);
title(ax,'Close-up around periapsis');

ax = nexttile(layout,3);
hold(ax,'on');
time = [min(vacuum.Time_s(1),drag.Time_s(1)),max(vacuum.Time_s(end),drag.Time_s(end))];
patch(ax,time([1 2 2 1]),band([1 1 2 2]),bandColor,'EdgeColor','none', ...
    'HandleVisibility','off');
plot(ax,drag.Time_s,vecnorm(drag.State(:,1:3),2,2)-radius,'Color',dragColor,'LineWidth',2.6);
plot(ax,vacuum.Time_s,vecnorm(vacuum.State(:,1:3),2,2)-radius, ...
    '--','Color',vacuumColor,'LineWidth',1.5);
yline(ax,band(2),'-','GRAM interface','LabelHorizontalAlignment','left');
yline(ax,band(1),':','Altitude floor','LabelHorizontalAlignment','left');
xlim(ax,time);
ylim(ax,[max(0,band(1)-0.08*diff(band)),band(2)+0.10*diff(band)]);
title(ax,'Altitude through the atmosphere');
timeLabel(ax,'Altitude (km)');

ax = nexttile(layout,5);
hold(ax,'on');
distance = vecnorm(drag.State(:,1:3),2,2);
entrySpeed = norm(vacuum.State(1,4:6));
entryDistance = norm(vacuum.State(1,1:3));
referenceSpeed = sqrt(entrySpeed^2+2*muVenus*(1./distance-1/entryDistance));
dragSpeed = vecnorm(drag.State(:,4:6),2,2);
yline(ax,0,'--','Color',vacuumColor,'DisplayName','Vacuum, no burn');
plot(ax,drag.Time_s,1000*(dragSpeed-referenceSpeed),'Color',dragColor, ...
    'LineWidth',2,'DisplayName','Drag on, no burn');
xlim(ax,time);
title(ax,'Speed loss at the same altitude');
timeLabel(ax,'Venus-relative speed difference (m/s)');

ax = nexttile(layout,6);
hold(ax,'on');
plot(ax,drag.Time_s,drag.DynamicPressure_Pa,'Color',dragColor,'LineWidth',2);
[peak,idx] = max(drag.DynamicPressure_Pa);
plot(ax,drag.Time_s(idx),peak,'o','MarkerFaceColor',dragColor,'MarkerEdgeColor','white');
title(ax,sprintf('Dynamic pressure  |  Peak %.3g Pa',peak));
xlim(ax,time);
timeLabel(ax,'Dynamic pressure (Pa)');

figures(2) = figure('Name','Earth–Venus–Earth post-flyby correction','Color','white', ...
    'Position',[90 110 1200 560]);
layout = tiledlayout(figures(2),1,2,'TileSpacing','compact','Padding','compact');
dates = string(datetime(mission.DatesJD,'ConvertFrom','juliandate','TimeZone','UTC','Format','yyyy-MM-dd'));
title(layout,strjoin(dates,'  →  '));

ax = nexttile(layout);
hold(ax,'on');
[~,~,au] = body_constants('Sun');
out = mission.Outbound.State(:,1:2)/au;
back = mission.VacuumReturn.State(:,1:2)/au;
plot(ax,out(:,1),out(:,2),'Color',[0.13 0.46 0.51],'LineWidth',2.2,'DisplayName','Earth to Venus');
returnPatch = mission.DragPatch;
if returnPatch.Status=="targeted"
    ret = returnPatch.ReturnState(:,1:2)/au;
    plot(ax,ret(:,1),ret(:,2),'Color',dragColor,'LineWidth',2.6, ...
        'DisplayName','Earth return after post-flyby burn');
end
plot(ax,back(:,1),back(:,2),'--','Color',vacuumColor,'LineWidth',1.5, ...
    'DisplayName','Nominal vacuum return');
plot(ax,0,0,'o','MarkerFaceColor',[1 0.75 0.29],'MarkerEdgeColor','none','DisplayName','Sun');
plot(ax,out(1,1),out(1,2),'o',back(end,1),back(end,2),'o', ...
    'Color',vacuumColor,'MarkerFaceColor','white','HandleVisibility','off');
plot(ax,out(end,1),out(end,2),'s','Color',dragColor,'MarkerFaceColor','white', ...
    'DisplayName','Venus encounter');
axis(ax,'equal');
grid(ax,'on');
title(ax,'Heliocentric transfer');
xlabel(ax,'Ecliptic X (AU)');
ylabel(ax,'Ecliptic Y (AU)');
legend(ax,'Location','southoutside','NumColumns',2);

ax = nexttile(layout);
hold(ax,'on');
drawVenus(ax,radius/1000,circle,'Venus');
vacuumPatch = mission.VacuumPatch;
if ~isempty(returnPatch.CoastState)
    coast = returnPatch.CoastState(:,1:3)*plane/1000;
    plot(ax,coast(:,2),coast(:,1),'Color',dragColor,'LineWidth',2.4, ...
        'DisplayName','Outbound coast after drag');
end
if ~isempty(vacuumPatch.CoastState)
    coastVac = vacuumPatch.CoastState(:,1:3)*plane/1000;
    plot(ax,coastVac(:,2),coastVac(:,1),'--','Color',vacuumColor,'LineWidth',1.5, ...
        'DisplayName','Outbound coast in vacuum');
end
if returnPatch.Status=="targeted"
    plot(ax,coast(end,2),coast(end,1),'d','MarkerFaceColor',dragColor, ...
        'MarkerEdgeColor','white','MarkerSize',9,'DisplayName','Post-flyby burn');
    subtitle(ax,sprintf('Post-flyby burn %.2f hours after nominal periapsis',returnPatch.CoastTime_s(end)/3600));
end
axis(ax,'equal');
grid(ax,'on');
title(ax,'Venus-centered coast to the burn');
xlabel(ax,'Along-track coordinate (10³ km)');
ylabel(ax,'Periapsis-axis coordinate (10³ km)');
legend(ax,'Location','southoutside','NumColumns',2);

    function drawFlyby(targetAxes)
        hold(targetAxes,'on');
        drawVenus(targetAxes,radius,circle,'Venus');
        n = size(circle,1);
        edge = (1:n-1).';
        ring = [(radius+band(2))*circle;(radius+band(1))*circle];
        patch(targetAxes,'Faces',[edge edge+1 n+edge+1 n+edge], ...
            'Vertices',ring,'FaceColor',bandColor,'EdgeColor','none','DisplayName',bandLabel);
        plot(targetAxes,approach(:,2),approach(:,1),'Color',[0.40 0.46 0.50], ...
            'LineWidth',1.7,'DisplayName','Shared inbound coast');
        plot(targetAxes,dragXY(:,2),dragXY(:,1),'Color',dragColor,'LineWidth',2.4, ...
            'DisplayName','Drag on, no burn');
        plot(targetAxes,vacXY(:,2),vacXY(:,1),'--','Color',vacuumColor,'LineWidth',1.35, ...
            'DisplayName','Vacuum, no burn');
        axis(targetAxes,'equal');
        grid(targetAxes,'on');
        xlabel(targetAxes,'Along-track coordinate (km)');
        ylabel(targetAxes,'Periapsis-axis coordinate (km)');
    end
end

function drawVenus(ax,radius,circle,label)
fill(ax,radius*circle(:,1),radius*circle(:,2),[0.93 0.87 0.76], ...
    'EdgeColor',[0.63 0.52 0.34],'LineWidth',0.8,'DisplayName',label,'Tag','VenusDisk');
end

function timeLabel(ax,verticalLabel)
xlabel(ax,'Seconds from nominal vacuum periapsis');
ylabel(ax,verticalLabel);
grid(ax,'on');
end
