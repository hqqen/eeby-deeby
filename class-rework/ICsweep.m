clear all; warning('off'); warning('this is an error with no valid ID')
%make sim objects
simWorld = World;
simWorld.groundPointsX = -2:12; simWorld.groundPointsY = -0.1*simWorld.groundPointsX;
simWorld.makeWorld();

biped = Robot;
% biped.q = [pi/96; 0; pi/24; -pi]; biped.dq = [1; 1.6; 1; -1.6];
biped.world = simWorld;
biped.B = [0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 1];
% q0 = [pi/12; 0; pi/12; 3.5*pi/4]; dq0 = [0;-pi/4;0;-pi/4];%[1; -1.6; 1; -1.6];

%def sim parameters
nSteps = 10; yS = [zeros(8,1)].'; yeS = []; tS = [0]; impacts = []; y = yS.';
opts = odeset("Events",@biped.isImpact);
tspan = [0 10];
rockTime = .1;

%create the coordinate space of ICs to iterate over (kick the swing leg to start)
dq = [0 0 0 .1].';
q0Space = [.01, pi/20:pi/20:pi]; %20 points in [0,pi] (even π curvature is probably too high )
tht0Space = 0:pi/40:pi/2;  %20 points in [0,pi/2] (anything greater than pi/2 corresponds to one leg starting up, less than zero is a leaning start)
q1Space = [.01, pi/20:pi/20:pi]; %dont let q0 or q1 be zero, singularities
tht1Space = 0:pi/40:pi/2;
qSpace = cell([max(size(q0Space)),max(size(tht0Space)),max(size(q1Space)),max(size(tht1Space))]);
numPoints = max(size(q0Space))*max(size(tht0Space))*max(size(q1Space))*max(size(tht1Space));
currPoint = 0;
%iterate through coordinate space
for q0Ind = 1:max(size(q0Space)) %doing it like this for a reason
    for tht0Ind = 1:max(size(tht0Space)) %doing it like this for a reason
        for q1Ind = 1:max(size(q1Space)) %doing it like this for a reason
            for tht1Ind = 1:max(size(tht1Space)) %doing it like this for a reason
                convergeFlag = 0; currPoint = currPoint + 1;
                %run the simulation for our chosen IC
                y = [q0Space(q0Ind); tht0Space(tht0Ind); q1Space(q1Ind); tht1Space(tht1Ind); dq];
                for i = 1:nSteps
                    clear te ye ie
%                     fprintf("Calculating Step " + i +"/" + nSteps + "\n");
                    [tTemp,yTemp,te,ye,ie] = ode45(@biped.doCtrl,tspan,y,opts);
%                     fprintf("Initial Foot Height:" + biped.getHeight(y(1:4)) + "\n")
                    yS = [yS;yTemp]; tS = [tS; tTemp + tS(end)]; impacts = [impacts, te]; yeS = [yeS; ye];
%                     fprintf("Step Time: " + tTemp(end) + "\n")
%                     fprintf("Final foot Height:" + biped.getHeight(yTemp(end,1:4)) + "\n")
                    

                    %check current warning state to see if there was a
                    %failure in the integrator; anything other than a 0x0
                    %char array is a failure
                    [~,ID] = lastwarn;
                    if (~isempty(ye) && isequal(ID,''))  %dont do impact dynamics if there is no impact
                        y = biped.doImpact(ye(1:4),ye(5:8));
                    else
%                         fprintf("Failed to find footfall event, ending integration... \n")
                        convergeFlag = 1;
                        break
                    end
                    %reset warning state
                    warning('this is an error with no valid ID that will not be raised to command line')
                end
                %assuming the simulation finished without running into
                %convergence issues, record stats at the corresponding
                %point in qSpace; if it did fail to converge write all NaNs
                %for processing later
                if convergeFlag == 0
                    %if the model succesfully finished integration record
                    %the average impact time, its standard deviation, the
                    %number of rocking footfalls (taking less than .05 sec)
                    %and the vector of impact times
                    qSpace(q0Ind,tht0Ind,q1Ind,tht1Ind) = {[mean(impacts);std(impacts);sum(impacts<.05,'all');sum(impacts < rockTime)]};%impacts};
                else
                    %else record a vector of NaNs that can be filtered out
                    %later (this makes indexing later less aggravating)
                    qSpace(q0Ind,tht0Ind,q1Ind,tht1Ind) = {[NaN;NaN;NaN;NaN;NaN]};
                end
                %progress indicator
                if mod(currPoint,100) == 0
                    fprintf("On Point " + currPoint + "/" + numPoints + " (" + round((currPoint/numPoints)*100,2) + " percent)\n")
                end
            end
        end
    end
end
