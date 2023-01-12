%y = [pi/12 -pi/6 pi/12 7*pi/6 0 0 0 0].';
y = [-pi/16 pi/12 -pi/16 5*pi/6 0 0 0 0].';
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
    fprintf("Step Time: " + tTemp(end) + "\n")
    if ~isempty(ye)  %dont do impact dynamics if there is no impact
        y(1) = ye(3);
        y(2) = ye(4) + pi; %w/o the +pi this works wonderfully
        y(3) = ye(1);
        y(4) = ye(2) + pi;
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

figure()
subplot(1,2,1); hold on;
plot(yS(:,2),yS(:,4)); 
xlabel("\theta_{st}"); ylabel("d\theta_{st}/dt")
subplot(1,2,2); hold on;
plot(yS(:,4),yS(:,8));
xlabel("\theta_{sw}"); ylabel("d\theta_{sw}/dt")
sgtitle("q_{st} =  q_{sw} = -\pi/16")

% for i = 1:max(size(yS))
% 
%     drawSelf(yS(i,:));
% %     exportgraphics(gcf,'gait.gif','Append',true);
%     % pause(.025);
% 
% end