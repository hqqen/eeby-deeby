%Gradient Descent Visualizer
%% User Defines
%Define the function to be zeroed
f = @(x,y) .1*x.^2 + .007*y.^4 + .004*y.^2 + .00005*x.^3 - .00234*x.^5 + .001*y.^6;
%Define the grid to work over - bounds defne where you think the zero is
%and stepsize effectives clarity of the contours (does not effect math tho)
gridSize = -10:.1:10;
%Set descent rate - higher is faster convergence in ideal cases but also
%more diverngences
gamma = .05;
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
contourf(gridX,gridY,f(real(gridX),real(gridY)),'ShowText',true); hold on;
contourf(gridX,gridY,f(real(gridX),real(gridY)), [eps,eps], 'r-', 'LineWidth',4, 'FaceColor', 'none');
contourf(gridX,gridY,f(real(gridX),real(gridY)), [-eps,-eps], 'r-', 'LineWidth',4, 'FaceColor', 'none');

%convert the provided function to symbolics, find its jacobian for use in
%calculating gradients
grad = jacobian(str2sym(char(f)));
delF = matlabFunction(grad);

%setup loop - make variables to track solver state
nStep = 0; nReset = 0;
while abs(f(x0(1),x0(2))) > eps

    %calculate and run update equation x(n+1) = x(n) - gamma*grad(f(x(n+1))
    x0Last = x0;
    x0 = x0 - gamma*g(x0(1),x0(2)).';

    if (x0(1) > gridSize(1) && x0(1) < gridSize(end)) && (x0(2) > gridSize(1) && x0(2) < gridSize(end))
        %if the new point is in the workspace put a red marker and draw a red
        %line from the last point to it
        plot(x0(1),x0(2),'r.',"MarkerSize",4);
        plot([x0Last(1), x0(1)],[x0Last(2), x0(2)],'r-')
    else
        %otherwise reinitialize the point in the workspace to avoid capture
        %by high gradients
        x0 = [(gridSize(end) - gridSize(1))*rand() + gridSize(1);...
            (gridSize(end) - gridSize(1))*rand() + gridSize(1)];
        plot(x0(1),x0(2),'g*',"MarkerSize",4)
        plot([x0Last(1), x0(1)],[x0Last(2), x0(2)],'g-')
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