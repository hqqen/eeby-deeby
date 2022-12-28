%Gradient Descent Visualizer
%% User Defines
%Define the function to be zeroed
f = @(x,y) .1*x.^2 + .007*y.^4 + .004*y.^2 + .00005*x.^3 - .00234*x.^5 + .001*y.^6;
%Define the grid to work over - bounds defne where you think the zero is
%and stepsize effectives clarity of the contours (does not effect math tho)
gridSize = -20:.1:20;
%Set descent rate - higher is faster convergence in ideal cases but also
%more diverngences
gamma = .1;
%Set solution tolerance (how close to zero the function needs to be to
%count as zeroed
eps = .1;

%% Sript Logic - Do not change!
%Randomly intialize guess in the workspace
x0 = [(gridSize(end) - gridSize(1))*rand() + gridSize(1);...
    (gridSize(end) - gridSize(1))*rand() + gridSize(1)];

%build contour plot of function; the target set is anywhere inside the
%thick red lnes
figure();
[gridX, gridY] = meshgrid(gridSize,gridSize);
contour3(gridX,gridY,f(real(gridX),real(gridY)),'ShowText',true); hold on;
contour3(gridX,gridY,f(real(gridX),real(gridY)), [eps,eps], 'p-', 'LineWidth',1, 'FaceColor', 'none');
contour3(gridX,gridY,f(real(gridX),real(gridY)), [-eps,-eps], 'p-', 'LineWidth',1, 'FaceColor', 'none');

%convert the provided function to symbolics, find its jacobian for use in
%calculating gradients
grad = jacobian(str2sym(char(f)));
delF = matlabFunction(grad);
g = matlabFunction(str2sym(char(f)));

%setup loop - make variables to track solver state
nStep = 0; nReset = 0;
xMax = [0,0,0]; xMin = [0,0,0];
while abs(f(x0(1),x0(2))) > eps

    %calculate and run update equation x(n+1) = x(n) - gamma*grad(f(x(n+1))
    x0Last = x0;
    x0 = x0 - gamma*delF(x0(1),x0(2)).';

%     %keep track of max and min values found in this solver step (theres
%     %probably a more elegant way to do this but w/e)
%     if x0(1) > xMax(1) && (x0(1) > gridSize(1) && x0(1) < gridSize(end))
%         xMax(1) = x0(1);
%     end
%     if x0(1) < xMin(1) && (x0(1) > gridSize(1) && x0(1) < gridSize(end))
%         xMin(1) = x0(1);
%     end
%     if x0(2) > xMax(2) && (x0(2) > gridSize(1) && x0(2) < gridSize(end))
%         xMax(2) = x0(2);
%     end
%     if x0(2) < xMin(2) && (x0(2) > gridSize(1) && x0(2) < gridSize(end))
%         xMin(2) = x0(2);
%     end
%     if g(x0(1),x0(2)) > xMax(3) && (g(x0(1),x0(2)) > gridSize(1) && g(x0(1),x0(2)) < gridSize(end))
%         xMax(3) = g(x0(1),x0(2));
%     end
%     if g(x0(1),x0(2)) < xMin(3) && (g(x0(1),x0(2)) > gridSize(1) && g(x0(1),x0(2)) < gridSize(end))
%         xMin(3) = g(x0(1),x0(2));
%     end
% 
%     xlim([xMin(1) - 10*eps, xMax(1) + 10*eps]);
%     ylim([xMin(2) - 10*eps, xMax(2) + 10*eps]);
%     zlim([xMin(3) - 10*eps, xMax(3) + 10*eps]);
    if (x0(1) > gridSize(1) && x0(1) < gridSize(end)) && (x0(2) > gridSize(1) && x0(2) < gridSize(end))
        %if the new point is in the workspace put a red marker and draw a red
        %line from the last point to it
        plot3(x0(1),x0(2),g(x0(1),x0(2)),'r.',"MarkerSize",4);
        plot3([x0Last(1), x0(1)],[x0Last(2), x0(2)],[g(x0Last(1),x0Last(2)), g(x0(1),x0(2))],'r-','LineWidth',2)
    else
        %otherwise reinitialize the point in the workspace to avoid capture
        %by high gradients
        x0 = [(gridSize(end) - gridSize(1))*rand() + gridSize(1);...
            (gridSize(end) - gridSize(1))*rand() + gridSize(1)];
        plot3(x0(1),x0(2),g(x0(1),x0(2)),'g*',"MarkerSize",4)
        plot3([x0Last(1), x0(1)],[x0Last(2), x0(2)],[g(x0Last(1),x0Last(2)), g(x0(1),x0(2))],'g-','LineWidth',2)
        nReset = nReset + 1;
    end
    pause(.2)

    %check if we've found a zero and print final outputs if we have, the
    %while loop breaks itself
    if abs(f(x0(1),x0(2))) < eps
        disp(nStep)
        disp(nReset)
        f(x0(1),x0(2))
    end
    nStep = nStep+1;
end