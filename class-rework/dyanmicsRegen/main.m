y = [pi/24 pi/6 pi/24 7*pi/6 0 0 0 0].';
tspan = [0 10];
nSteps = 10;
opts = odeset("Events",@isImpact);

yS = []; tS = [0]; yeS = [];

for i = 1:nSteps

    fprintf("Calculating Step " + i +"/" + nSteps + "\n");
    [tTemp,yTemp,te,ye,~] = ode45(@fdbkCtrl,tspan,y,opts);
    yS = [yS;yTemp]; 
    tS = [tS; tTemp + tS(end)]; 
    yeS = [yeS; ye];
    %     currStep = [currStep; ones(size(tTemp))*i]; %so plotted knows what step we're on
    fprintf("Step Time: " + tTemp(end) + "\n")
    if ~isempty(ye)  %dont do impact dynamics if there is no impact
        y(1) = ye(3);
        y(2) = ye(4) - ye(1) + pi; %w/o the +pi this works wonderfully
        y(3) = ye(1);
        y(4) = ye(2) + ye(1) + pi;
        y(5) = ye(7);
        y(6) = ye(8);
        y(7) = ye(5);
        y(8) = ye(6);
        fprintf("Step " + i + " Complete! \n")
    else
        fprintf("Failed to find footfall event, ending integration... \n")
        break
    end

end