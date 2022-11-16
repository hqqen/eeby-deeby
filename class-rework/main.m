clear all;
%make sim objects
simWorld = World;
simWorld.groundPointsX = -2:12; simWorld.groundPointsY = -0.1*simWorld.groundPointsX;
simWorld.makeWorld();

biped = Robot;
% biped.q = [pi/96; 0; pi/24; -pi]; biped.dq = [1; 1.6; 1; -1.6];
biped.world = simWorld;
biped.B = [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1];
q0 = [pi/12; pi/12; pi/4; pi/2]; dq0 = [0;0;0;.1];%[1; -1.6; 1; -1.6];

%def sim parameters
nSteps = 10; yS = [q0;dq0].'; yeS = []; tS = [0]; impacts = [0]; y = yS.';
opts = odeset("Events",@biped.isImpact);
tspan = [0 10];
rockTime = .25; %any step shorter than this is a rocking motion
currStep = [];

%run sim
for i = 1:nSteps
    fprintf("Calculating Step " + i +"/" + nSteps + "\n");
    [tTemp,yTemp,te,ye,~] = ode45(@biped.doCtrl,tspan,y,opts);
    fprintf("Initial Foot Height:" + biped.getHeight(y(1:4)) + "\n")
    yS = [yS;yTemp]; tS = [tS; tTemp + tS(end)]; impacts = [impacts, te(end)]; yeS = [yeS; ye];
%     currStep = [currStep; ones(size(tTemp))*i]; %so plotted knows what step we're on
    fprintf("Step Time: " + tTemp(end) + "\n")
    fprintf("Final foot Height:" + biped.getHeight(yTemp(end,1:4)) + "\n")
    if ~isempty(ye)  %dont do impact dynamics if there is no impact
        y = biped.doImpact(ye(1:4),ye(5:8));
    else
        fprintf("Failed to find footfall event, ending integration... \n")
        break
    end
end
if any(impacts < rockTime)
    fprintf("Rocking Motion Detected \n")
    fprintf(sum(impacts < rockTime) + " Rocking Steps Detected \n")
end
biped.drawStates(yS,tS,ye,te)
biped.drawHeight(yS,tS)
%RUN ON 2022a OR LATER!!!!!!!!!!!!!
%exportgraphics(biped.animate(yS,tS),'firstSteps.gif')